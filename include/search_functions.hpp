/*
 * search_functions.hpp
 *
 *  Created on: 2018年12月22日
 *      Author: sparklexs
 */

#ifndef INCLUDE_SEARCH_FUNCTIONS_HPP_
#define INCLUDE_SEARCH_FUNCTIONS_HPP_

#include "k-ary_search_tree.hpp"
#include <immintrin.h>

namespace karytree {
/*
 * FIXME: a really bad implementation, works totally astray
 * @deprecated: biased or skewed search
 * @endindex = size - 1
 */
template<indextype k>
indextype search_sorted_array(keytype goal, const keytype *target,
		indextype endindex) {
	if (_UNLIKELY(endindex <= 1)) {
		//1 or 2 elements
		if (*target < goal)
			return endindex;
		else
			return 0;
	}
	// now we have height >=2
	indextype height = ktree<k>::getTreeHeight(endindex + 1);
	indextype remainder = (endindex + 1) - ktree<k>::PERFECT_SIZES[height];
	indextype depth = 0, left = 0;

	if (_LIKELY(remainder)) {
		remainder = (endindex + 1) - ktree<k>::PERFECT_SIZES[height - 1];
		// the following search seems to be an uneven, or biased search
		// it first search the front @remainder elements, and then
		// the remaining ones.
		indextype step = -1, offset;
		for (uint32_t i = 1; i < k; i++) {
			offset = i * (remainder - 1) / (k - 1);
			if (goal == target[offset])
				return offset;
			else if (goal < target[offset])
				break;
			step = offset;
		}
		left = step + 1;
		depth++;
	}

	for (; depth < height; ++depth) {
		indextype step = -1, offset;
		for (indextype i = 1; i < k; ++i) {
			// here "+1" aovids PERFECT_SIZES[0] == 0
			offset = i * (ktree<k>::PERFECT_SIZES[height - depth - 1] + 1) - 1;
			indextype separatrorIndex = left + offset;
			// we can combine these two comparison into one
			// and add one comparison with left+1 before return
			// but seems inefficient
			if (goal == target[separatrorIndex])
				return separatrorIndex;
			else if (goal < target[separatrorIndex])
				break;
			step = offset;
		}
		left += step + 1;
	}
	return left;
}

template<indextype k>
indextype search_sorted_array2(keytype goal, const keytype *target,
		indextype endindex) {
	if (_UNLIKELY(endindex <= 1)) {
		//1 or 2 elements
		if (*target < goal)
			return endindex;
		else
			return 0;
	}
	// now we have height >=2
	indextype height = ktree<k>::getTreeHeight(endindex + 1);
	indextype remainder = (endindex + 1) - ktree<k>::PERFECT_SIZES[height];
	indextype depth = 0, left = 0;

	if (_LIKELY(remainder)) {
		indextype interval_edge = ktree<k>::PERFECT_SIZES[height - 1];
		indextype interval_mid = (endindex + 1
				- 2 * ktree<k>::PERFECT_SIZES[height - 1] + k - 3) / (k - 2);

		indextype os = interval_edge;
		if (goal == target[os])
			return os;
		else if (goal < target[os]) {
			left = 0;
			goto nextlevel;
		}

		for (uint32_t i = 0; i < k - 2; i++) {
			os += interval_mid;
			if (goal == target[os])
				return os;
			else if (goal < target[os]) {
				left = os - interval_mid;
				goto nextlevel;
			}
		}
		left = endindex - interval_edge;
		nextlevel: depth++;
	}

	for (; depth < height; ++depth) {
		indextype step = -1, offset;
		for (indextype i = 1; i < k; ++i) {
			// here "+1" aovids PERFECT_SIZES[0] == 0
			offset = i * (ktree<k>::PERFECT_SIZES[height - depth - 1] + 1) - 1;
			indextype separatrorIndex = left + offset;
			// we can combine these two comparison into one
			// and add one comparison with left+1 before return
			// but seems inefficient
			if (goal == target[separatrorIndex])
				return separatrorIndex;
			else if (goal < target[separatrorIndex])
				break;
			step = offset;
		}
		left += step + 1;
	}
	return left;

}

// FIXME: assume each node is full of k-1 elements
// this function should be moved into the ktree class
template<indextype k>
indextype search_linearized_tree(int *foundp, keytype goal,
		const keytype *target, indextype endindex) {
	indextype left = 0, next = 0;
	indextype m;
	while (next < endindex) {
		left = next;
		for (m = 0; m < k - 1; m++) {
			if (target[left + m] == goal) {
				*foundp = 1;
				return left + m;
			}
			if (target[left + m] > goal)
				break;
		}
		next = (left + 1) * k + m * (k - 1) - 1;
	}
	// now we have left > goal and m in range [0,k-1]
	*foundp = 0;
	// m might be k -1, which means the largest value of the tree is less than goal
	return left + m;
}

template<>
indextype search_linearized_tree<3>(int *foundp, keytype goal,
		const keytype *target, indextype endindex) {
	indextype left = 0, next = 0;
	indextype m = 0;
	while (next < endindex) {
		left = next;
		__m64 pgoal = _mm_set1_pi32(goal);
		__m64 ptarget = _mm_setr_pi32(target[left], target[left + 1]);
		__m64 presult = _mm_cmpeq_pi32(pgoal, ptarget);
		int result = _mm_movemask_pi8(presult);
		if (result) {
			*foundp = 1;
			return left + __builtin_ctz(result);
		}
		presult = _mm_cmpgt_pi32(ptarget, pgoal);
		result = _mm_movemask_pi8(presult);
		if (result)
			m = __builtin_ctz(result) >> 2;
		else
			m = 2;

		next = left * 3 + (m << 1) + 2;
	}
	// now we have left > goal and m in range [0,k-1]
	*foundp = 0;
	// m might be k -1, which means the largest value of the tree is less than goal
	return left + m;
}

template<>
indextype search_linearized_tree<5>(int *foundp, keytype goal,
		const keytype *target, indextype endindex) {
	indextype left = 0, next = 0;
	indextype m = 0;
	while (next < endindex) {
		left = next;
		__m128i pgoal = _mm_set1_epi32(goal);
		__m128i ptarget = _mm_lddqu_si128((__m128i *) target);
		__m128i presult = _mm_cmpeq_epi32(pgoal, ptarget);
		int result = _mm_movemask_ps(_mm_castsi128_ps(presult));
		if (result) {
			*foundp = 1;
			return left + __builtin_ctz(result);
		}
		presult = _mm_cmpgt_epi32(ptarget, pgoal);
		result = _mm_movemask_ps(_mm_castsi128_ps(presult));
		if (result)
			m = __builtin_ctz(result);
		else
			m = 4;

		next = left * 5 + (m << 2) + 4;
	}
	// now we have left > goal and m in range [0,k-1]
	*foundp = 0;
	// m might be k -1, which means the largest value of the tree is less than goal
	return left + m;
}

template<>
indextype search_linearized_tree<9>(int *foundp, keytype goal,
		const keytype *target, indextype endindex) {
	indextype left = 0, next = 0;
	indextype m = 0;
	while (next < endindex) {
		left = next;
		__m256i pgoal = _mm256_set1_epi32(goal);
		__m256i ptarget = _mm256_lddqu_si256((__m256i *) target);
		__m256i presult = _mm256_cmpeq_epi32(pgoal, ptarget);
		int result = _mm256_movemask_ps(_mm256_castsi256_ps(presult));
		if (result) {
			*foundp = 1;
			return left + __builtin_ctz(result);
		}
		presult = _mm256_cmpgt_epi32(ptarget, pgoal);
		result = _mm256_movemask_ps(_mm256_castsi256_ps(presult));
		if (result)
			m = __builtin_ctz(result);
		else
			m = 8;

		next = left * 9 + (m << 3) + 8;
	}
	// now we have left > goal and m in range [0,k-1]
	*foundp = 0;
	// m might be k -1, which means the largest value of the tree is less than goal
	return left + m;
}

// FIXME: assume each node is full of k-1 elements
template<indextype k>
indextype search_linearized_tree_bounded(int *foundp, keytype goal,
		const keytype *target, indextype &depth, indextype &rank,
		indextype endindex) {
	indextype next = ktree<k>::PERFECT_SIZES[depth] + rank * (k - 1), left;
	indextype m;
	do {
		left = next;
		for (m = 0; m < k - 1; m++) {
			if (target[left + m] == goal) {
				*foundp = 1;
				return left + m;
			}
			if (target[left + m] > goal)
				break;
		}
		depth++;
		rank = rank * k + m;
		// next = (left + 1) * k + m * (k - 1) - 1;
		next = ktree<k>::PERFECT_SIZES[depth] + rank * (k - 1);
	} while (next < endindex);

	depth--;
	rank = (rank - m) / k;

	*foundp = 0;
	return left + m;
}

// FIXME: assume each node is full of k-1 elements
template<>
indextype search_linearized_tree_bounded<3>(int *foundp, keytype goal,
		const keytype *target, indextype &depth, indextype &rank,
		indextype endindex) {
	indextype next = ktree<3>::PERFECT_SIZES[depth] + (rank << 1), left;
	indextype m;
	do {
		left = next;

		__m64 pgoal = _mm_set1_pi32(goal);
		__m64 ptarget = _mm_setr_pi32(target[left], target[left + 1]);
		__m64 presult = _mm_cmpeq_pi32(pgoal, ptarget);
		int result = _mm_movemask_pi8(presult);
		if (result) {
			*foundp = 1;
			return left + (result >> 7);
		}
		presult = _mm_cmpgt_pi32(ptarget, pgoal);
		result = _mm_movemask_pi8(presult);
		if (result)
			m = __builtin_ctz(result) >> 2;
		else
			m = 2;

		depth++;
		rank = rank * 3 + m;
		// next = (left + 1) * k + m * (k - 1) - 1;
		next = ktree<3>::PERFECT_SIZES[depth] + (rank << 1);
	} while (next < endindex);

	depth--;
	rank = (rank - m) / 3;

	*foundp = 0;
	return left + m;
}

// FIXME: assume each node is full of k-1 elements
template<>
indextype search_linearized_tree_bounded<5>(int *foundp, keytype goal,
		const keytype *target, indextype &depth, indextype &rank,
		indextype endindex) {
	indextype next = ktree<5>::PERFECT_SIZES[depth] + (rank << 2), left;
	indextype m;
	do {
		left = next;

		__m128i pgoal = _mm_set1_epi32(goal);
		__m128i ptarget = _mm_lddqu_si128((__m128i *) target);
		__m128i presult = _mm_cmpeq_epi32(pgoal, ptarget);
		int result = _mm_movemask_ps(_mm_castsi128_ps(presult));
		if (result) {
			*foundp = 1;
			return left + __builtin_ctz(result);
		}
		presult = _mm_cmpgt_epi32(ptarget, pgoal);
		result = _mm_movemask_ps(_mm_castsi128_ps(presult));
		if (result)
			m = __builtin_ctz(result);
		else
			m = 4;

		depth++;
		rank = rank * 5 + m;
		// next = (left + 1) * k + m * (k - 1) - 1;
		next = ktree<5>::PERFECT_SIZES[depth] + (rank << 2);
	} while (next < endindex);

	depth--;
	rank = (rank - m) / 5;

	*foundp = 0;
	return left + m;
}

// FIXME: assume each node is full of k-1 elements
template<>
indextype search_linearized_tree_bounded<9>(int *foundp, keytype goal,
		const keytype *target, indextype &depth, indextype &rank,
		indextype endindex) {
	indextype next = ktree<9>::PERFECT_SIZES[depth] + (rank << 3), left;
	indextype m;
	do {
		left = next;

		__m256i pgoal = _mm256_set1_epi32(goal);
		__m256i ptarget = _mm256_lddqu_si256((__m256i *) target);
		__m256i presult = _mm256_cmpeq_epi32(pgoal, ptarget);
		int result = _mm256_movemask_ps(_mm256_castsi256_ps(presult));
		if (result) {
			*foundp = 1;
			return left + __builtin_ctz(result);
		}
		presult = _mm256_cmpgt_epi32(ptarget, pgoal);
		result = _mm256_movemask_ps(_mm256_castsi256_ps(presult));
		if (result)
			m = __builtin_ctz(result);
		else
			m = 8;

		depth++;
		rank = rank * 9 + m;
		// next = (left + 1) * k + m * (k - 1) - 1;
		next = ktree<9>::PERFECT_SIZES[depth] + (rank << 3);
	} while (next < endindex);

	depth--;
	rank = (rank - m) / 9;

	*foundp = 0;
	return left + m;
}

}

#endif /* INCLUDE_SEARCH_FUNCTIONS_HPP_ */
