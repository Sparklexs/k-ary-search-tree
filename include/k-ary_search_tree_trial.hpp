/*
 * k-ary_search_tree_trial.cpp
 *
 *  Created on: 2018Äê1ÔÂ15ÈÕ
 *      Author: John
 */

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <iostream>

#include "synthetic.h"
#include "timer.h"

// specially works for k equals 3
namespace karytree_trial {
typedef uint32_t keytype;
typedef size_t indextype;

const static indextype PERFECT_SIZES[] = { 0, 2, 8, 26, 80, 242, 728, 2186,
		6560, 19682, 59048, 177146, 531440, 1594322, 4782968, 14348906,
		43046720, 129140162, 387420488, 1162261466, 3486784400 /*, 10460353202, 31381059608,
 94143178826, 282429536480, 847288609442, 2541865828328, 7625597484986,
 22876792454960, 68630377364882, 205891132094648, 617673396283946,
 1853020188851840, 5559060566555522*/};

const static indextype K_s[] = { 1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683,
		59049, 177147, 531441, 1594323, 4782969, 14348907, 43046721, 129140163,
		387420489, 1162261467, 3486784401, 10460353203, 31381059609,
		94143178827, 282429536481, 847288609443, 2541865828329, 7625597484987,
		22876792454961, 68630377364883, 205891132094649, 617673396283947,
		1853020188851841, };

// the static array PERFECT_SIZES is only used for k equals 3
// so we have to calculate actual array lengths when k varies.
template<uint32_t k = 3>
inline void initializeArrayLength(std::vector<keytype>& vec) {
	vec[0] = 0;
	vec[1] = k - 1;
	for (indextype i = 1; i < vec.size(); i++) {
		vec[i] = (vec[i - 1] + 1) * k - 1;
	}
}

/*
 * @return height are always larger than 1
 */
indextype getTreeHeight(indextype size) {
	indextype height = 0;
	while (PERFECT_SIZES[height] < size)
		height++;
	return height;
}

/////////////////////////////////////////////////////////////////

/*
 * @endindex = size - 1
 */
template<indextype k = 3>
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
	indextype height = getTreeHeight(endindex + 1);
	indextype remainder = (endindex + 1) - PERFECT_SIZES[height];
	indextype depth = 0, left = 0;

	if (_LIKELY(remainder)) {
		remainder = (endindex + 1) - PERFECT_SIZES[height - 1];
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
			offset = i * (PERFECT_SIZES[height - depth - 1] + 1) - 1;
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

// assume each node is full of k-1 elements
// this function should be moved into the ktree class
template<indextype k = 3>
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

template<indextype k = 3>
indextype search_linearized_tree_bounded(int *foundp, keytype goal,
		const keytype *target, indextype &depth, indextype &rank,
		indextype endindex) {
	indextype next = PERFECT_SIZES[depth] + rank * (k - 1), left;
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
		next = PERFECT_SIZES[depth] + rank * (k - 1);
	} while (next < endindex);

	depth--;
	rank = (rank - m) / k;

	*foundp = 0;
	return left + m;
}

// always start search from root and return target (depth,rank)
template<indextype k = 3>
indextype search_linearized_tree_with_result(int *foundp, keytype goal,
		const keytype *target, indextype &depth, indextype &rank,
		indextype endindex) {
	indextype next = 0, left;
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
		next = PERFECT_SIZES[depth] + rank * (k - 1);
	} while (next < endindex);

	depth--;
	rank = (rank - m) / k;

	// m might be k -1, which means the largest value of the tree is less than goal
	*foundp = 0;
	return left + m;
}

/////////////////////////////////////////////////////////////////

/*
 * deprecated!
 * @param index is used for sorted array
 * @param pos is used for linearized ktree
 */
//template<indextype k = 3>
//indextype index_to_pos(indextype index, indextype size) {
//	indextype height = getTreeHeight(size);
//	indextype depth = 0, depth_prev = 0, offset = 0;
//	while (true) {
//		indextype maxheight = 1;
//		while (index >= PERFECT_SIZES[maxheight]) {
//			maxheight++;
//		}
//		depth = height - maxheight;
//		// when gaps between depth and depth_prev is larger than 2,
//		// how to transfer the offset into current level?
//		// seems we have to store all the depth_prev to implement
//		// this goal
//		if (depth < depth_prev + 1)
//			offset += index / (PERFECT_SIZES[maxheight - 1] + 1);
//		else
//			offset = offset
//					* (PERFECT_SIZES[depth - depth_prev]
//							- PERFECT_SIZES[depth - depth_prev - 1])
//					+ index / (PERFECT_SIZES[maxheight - 1] + 1);
//		indextype remainder = index % (PERFECT_SIZES[maxheight - 1] + 1);
//		if (remainder == PERFECT_SIZES[maxheight - 1]) {
//			break;
//		}
//		index = remainder;
//		depth_prev = depth;
//	}
//}
/*
 * given an index in a sorted array, find its
 * corresponding position in the linearized search tree
 * @param index is used for sorted array
 * @param pos is used for linearized ktree
 */
template<indextype k = 3>
indextype index_to_pos(indextype index, indextype size) {
	indextype height = getTreeHeight(size);
	indextype remainder = (size) - PERFECT_SIZES[height];

	// complete tree
	if (_LIKELY(remainder)) {
		remainder = (size) - PERFECT_SIZES[height - 1];
		indextype fringe_entry = (remainder - 1) * k / (k - 1);
		if (index > fringe_entry) {
			index = index - remainder;
			height--;
		}
	}
	// perfect tree
	indextype digits[100];
	indextype msb = 0;
	indextype ind = index;
	while (ind != 0) {
		indextype rmd = ind % k;
		ind = ind / k;
		digits[msb] = rmd;
		msb++;
	}
	digits[msb] = 0;
	// pos--;
	indextype lsb = 0;
	while (digits[lsb] == k - 1)
		lsb++;
	indextype depth = height - lsb - 1;
	indextype offset = 0;
	for (indextype i = msb - 1; i > lsb; i--) {
		offset = offset * k + digits[i] * (k - 1);
	}
	offset += digits[lsb] + PERFECT_SIZES[depth];
	return offset;
}

/*
 * given a position in a linearized ktree, find its
 * corressponding index in the sorted array
 * @param pos is used for linearized ktree
 * @param index is used for sorted array
 */
template<indextype k = 3>
indextype pos_to_index(indextype pos, indextype size) {
	indextype height = getTreeHeight(size);

	indextype depth = 1;
	while (pos >= PERFECT_SIZES[depth]) {
		depth++;
	}
	depth--;
	indextype offset = pos - PERFECT_SIZES[depth];
	// three parts
	// 1. its left sibling trees (+1 means the elements at its same level)
	// 2. sub trees under it
	// 3. offset of itself inside the node composed of k-1 elements
	indextype index = offset / (k - 1) * (PERFECT_SIZES[height - depth] + 1)
			+ (offset % (k - 1) + 1) * PERFECT_SIZES[height - depth - 1]
			+ offset % (k - 1);
	if (index < size)
		return index;
	else {
		indextype remainder = (size) - PERFECT_SIZES[height - 1];
//		indextype fringe_entry = (remainder - 1) * k / (k - 1);
		height--;
		return offset / (k - 1) * (PERFECT_SIZES[height - depth] + 1)
				+ (offset % (k - 1) + 1) * PERFECT_SIZES[height - depth - 1]
				+ offset % (k - 1) + remainder;

	}
}

/*
 * given a position in the linearized ktree
 * find the index of its right-hand value
 * in the sorted array
 */
template<indextype k = 3>
indextype next_index(indextype pos, indextype size) {
	return pos_to_index<k>(pos, size) + 1;
}

/*
 * given a position in the linearized ktree
 * find the position of its right-hand value
 * from the sorted array
 */
template<indextype k = 3>
indextype next_pos(indextype pos, indextype size) {
	indextype height = getTreeHeight(size);
	indextype remainder = (size) - PERFECT_SIZES[height];
	indextype depth_next = 0, offset_next = 0;

	indextype depth = 1;
	while (pos >= PERFECT_SIZES[depth]) {
		depth++;
	}
	depth--;
	indextype offset = pos - PERFECT_SIZES[depth];

	// complete tree
	if (_LIKELY(remainder)) {
		// here "- 1" means the offset
		// remainder also means the offset of the last leaf node
		remainder = (size) - PERFECT_SIZES[height - 1] - 1;

		if (depth + 1 < height) {
			// elements that are not in leaves, their next
			// elements will certainly fall into leaf nodes
			offset_next = offset / (k - 1) * (K_s[height - depth - 1] * (k - 1))
					+ (offset % (k - 1) + 1) * K_s[height - depth - 2]
							* (k - 1);
			if (offset_next > remainder) {
				height--;
				goto PERFECT_TREE;
			}
			return offset_next + PERFECT_SIZES[height - 1];
		} else {
			// for those fall into leaf-level
			if (offset == remainder) {
				// handle this specially
				offset = (offset + k - 2) / (k - 1) * (k - 1) - 1;
				goto ASCEND;
			}
			goto LEAF;
		}
	} else {
		// perfect tree
		PERFECT_TREE: if (depth + 1 < height) {
			// elements that are not in leaves, their next
			// elements will certainly fall into leaf nodes
			depth_next = height - 1;
			offset_next = (offset / (k - 1)) * K_s[height - depth - 1] * (k - 1)
					+ (offset % (k - 1) + 1) * K_s[height - depth - 2]
							* (k - 1);
			return offset_next + PERFECT_SIZES[depth_next];
		}
		// only elements from leaf nodes, their next elements can be
		// right next to them or grow into upper nodes
		LEAF: if ((offset + 1) % (k - 1) != 0) {
			return pos + 1;
		} else if (offset == K_s[height - 1] * (k - 1) - 1) {
			// manually circular-join the last leaf to the root
			return 0;
		} else {
			ASCEND: indextype sum_of_middle_points = 0;
			while ((offset + 1) % (K_s[height - depth_next - 2] * (k - 1)) != 0) {
				sum_of_middle_points += (offset + 1)
						/ (K_s[height - depth_next - 2] * (k - 1));
				depth_next++;
			}
			offset_next = (offset + 1)
					/ (K_s[height - depth_next - 2] * (k - 1))
					- sum_of_middle_points - 1;
			return offset_next + PERFECT_SIZES[depth_next];
		}
	}
}

// sequentially output the values in ascending order
template<indextype k = 3>
void traverse_tree(keytype *pointer, indextype rank, indextype depth,
		indextype size) {
	for (indextype i = 0; i < k; i++) {
		if (PERFECT_SIZES[depth + 1] + (k * rank + i) * (k - 1) < size) {
			traverse_tree(pointer, rank * k + i, depth + 1, size);
		}
		if (i != k - 1)
			std::cout << pointer[PERFECT_SIZES[depth] + (k - 1) * rank + i]
					<< ",";
	}
	std::cout << std::endl;
}

/////////////////////////////////////////////////////////////////

template<indextype k = 3>
void construct_tree(keytype* array, indextype size) {
	indextype height = getTreeHeight(size);
	indextype remainder = (size) - PERFECT_SIZES[height];
	std::vector<keytype> tmp(size);
	indextype digits[100];
	tmp[PERFECT_SIZES[height - 1]] = array[0];

// complete tree
	if (_LIKELY(remainder)) {
		remainder = (size) - PERFECT_SIZES[height - 1];
		indextype fringe_entry = (remainder - 1) * k / (k - 1);

		for (indextype i = 1; i <= fringe_entry; i++) {
			indextype pos = 0;
			indextype ind = i;
			while (ind != 0) {
				indextype rmd = ind % k;
				ind = ind / k;
				digits[pos] = rmd;
				pos++;
			}
			digits[pos] = 0;
			// pos--;
			indextype j = 0;
			while (digits[j] == k - 1)
				j++;
			indextype depth = height - j - 1;
			indextype offset = 0;
			for (indextype p = pos - 1; p > j; p--) {
				offset = offset * k + digits[p] * (k - 1);
			}
			offset += digits[j] + PERFECT_SIZES[depth];
			tmp[offset] = array[i];
		}
		for (indextype i = fringe_entry + 1; i < size; i++) {
			indextype pos = 0;
			indextype ind = i - remainder;
			while (ind != 0) {
				indextype rmd = ind % k;
				ind = ind / k;
				digits[pos] = rmd;
				pos++;
			}
			digits[pos] = 0;
			// pos--;
			indextype j = 0;
			while (digits[j] == k - 1)
				j++;
			indextype depth = height - j - 2;
			indextype offset = 0;
			for (indextype p = pos - 1 > pos ? 0 : pos - 1; p > j; p--) {
				offset = offset * k + digits[p] * (k - 1);
			}
			offset += digits[j] + PERFECT_SIZES[depth];
			tmp[offset] = array[i];
		}
	}
// perfect tree
	else {
		for (indextype i = 1; i < size; i++) {
			indextype pos = 0;
			indextype ind = i;
			while (ind != 0) {
				indextype rmd = ind % k;
				ind = ind / k;
				digits[pos] = rmd;
				pos++;
			}
			// in case the previous ind has more digits
			digits[pos] = 0;
			// pos--;
			indextype j = 0;
			while (digits[j] == k - 1)
				j++;
			indextype depth = height - j - 1;
			indextype offset = 0;
			for (indextype p = pos - 1; p > j; p--) {
				offset = offset * k + digits[p] * (k - 1);
			}
			offset += digits[j] + PERFECT_SIZES[depth];
			tmp[offset] = array[i];
		}
	}
	memmove(array, tmp.data(), size * sizeof(keytype));
}

/*
 * for nodes in the perfect tree
 * @param index_array points to the start position of current subtree
 * @param index_tree points to the very head of the root of the whole tree
 */
template<indextype k = 3>
void reorder(keytype *index_array, indextype rank, indextype depth,
		indextype height, keytype *index_tree) {
	for (indextype i = 0; i < k - 1; i++) {
		index_tree[PERFECT_SIZES[depth] + (k - 1) * rank + i] =
				index_array[PERFECT_SIZES[height - 1] * (i + 1) + i];
	}
	if (--height) {
		// FIXME: 1. for the leaf leaf level
		// we have only k-1 elements left
		// and they can be sequentially
		// copied to destination address
		// 2. index_tree should be also
		// updated for less computation
		// 3. merge this for-loop with
		// the previous one?

		for (indextype i = 0; i < k; i++) {
			reorder<k>(index_array + i * (PERFECT_SIZES[height] + 1),
					rank * k + i, depth + 1, height, index_tree);
		}
	}
}

template<indextype k = 3>
void reorder_with_remainder(keytype *index_array, indextype size,
		indextype remainder, indextype rank, indextype depth, indextype height,
		keytype *index_tree) {
	// @branch denotes the branch where lies the complete subtree
	indextype branch = remainder / (K_s[height - 2] * (k - 1));
	//if no complete tree exists, there will have a overlapping of one parent node

	// 1. the left sibling nodes are perfect trees with height of h-1
	for (indextype i = 0; i < branch; i++) {
		// here subtrees go before the parent nodes
		reorder<k>(index_array + i * (PERFECT_SIZES[height - 1] + 1),
				rank * k + i/*rank*/, depth + 1/*depth*/, height - 1,
				index_tree);
		index_tree[PERFECT_SIZES[depth] + rank * (k - 1) + i] =
				index_array[PERFECT_SIZES[height - 1] * (i + 1) + i];
	}

	// 2. the 'branch'-th subtree is one complete tree with height of h-1
	// or h-2 and it only exists when following condition is true
	// FIXME: what if complete tree becomes a perfect tree of h-2?
	// the remainder becomes 0 and are the following codes compatible
	// with such situation?
	if (remainder % ((K_s[height - 1] * (k - 1)) / k)) {
		indextype next_remainder = remainder
				% ((K_s[height - 1] * (k - 1)) / k);
		if (height == 2) {
			for (indextype i = 0; i < next_remainder; i++) {
				index_tree[size - next_remainder + i] =
						index_array[PERFECT_SIZES[1] * branch + branch + i];
			}
		} else {
			reorder_with_remainder<k>(
					index_array + branch * (PERFECT_SIZES[height - 1] + 1),
					size, next_remainder, rank * k + branch, depth + 1,
					height - 1, index_tree);
		}
		branch++;
	}
	keytype* right_sibling_start = index_array
			+ PERFECT_SIZES[height - 2] * branch + branch - 1 + remainder;
	// 3. the right sibling subtrees are perfect trees with height of h-2
	branch--;
	for (indextype i = branch; i < k - 1; i++) {
		// here is a little different from the first part, as the parent nodes
		// go before the subtrees, so rank*k + i '+ 1'
		index_tree[PERFECT_SIZES[depth] + rank * (k - 1) + i] =
				right_sibling_start[(i - branch) * PERFECT_SIZES[height - 2]
						+ (i - branch)];
		if (height != 2)
			reorder<k>(
					right_sibling_start
							+ (i - branch) * PERFECT_SIZES[height - 2]
							+ (i + 1 - branch), rank * k + i + 1/*rank*/,
					depth + 1, height - 2, index_tree);
	}
}

template<indextype k = 3>
void construct_tree_fast(keytype* array, indextype size) {
	indextype height = getTreeHeight(size);
	indextype remainder = (size) - PERFECT_SIZES[height];
	std::vector<keytype> tmp(size);

// complete tree
	if (_LIKELY(remainder)) {
		remainder = (size) - PERFECT_SIZES[height - 1];
		// remainder also means the number of elements in
		// the leaf level
		reorder_with_remainder<k>(array, size, remainder, 0, 0, height,
				tmp.data());
	}
// perfect tree
	else {
		reorder<k>(array, 0, 0, height, tmp.data());
	}
	memmove(array, tmp.data(), size * sizeof(keytype));
}

/////////////////////////////////////////////////////////////////

template<indextype k = 3>
void intersect_sequential(keytype *small, indextype size_small, keytype *large,
		indextype size_large, std::vector<keytype> & out) {
	out.clear();
	// linearly fetch
	int found = 0;
	for (indextype i = 0; i < size_small; i++) {
		/*indextype pos = */search_linearized_tree<k>(&found, small[i], large,
				size_large - 1);
//		std::cout << i << ": " << small[i] << ", " << large[pos] << std::endl;
		if (found)
			out.emplace_back(small[i]);
	}
}

template<indextype k = 3>
void intersect_sorted(keytype *small, indextype size_small, keytype *large,
		indextype size_large, std::vector<keytype> & out) {
	out.clear();
	// fetch by value
	int found = 0;
	std::function<void(indextype rank, indextype depth)> search;
	search =
			[&]( indextype rank, indextype depth)-> void {
				for (indextype i = 0; i < k; i++) {
					if (PERFECT_SIZES[depth + 1] + (k * rank + i) * (k - 1) < size_small) {
						search(rank * k + i, depth + 1);
					}
					if (i != k - 1) {
						indextype pos = search_linearized_tree<k>(&found,
								small[PERFECT_SIZES[depth] + (k - 1) * rank + i],
								large, size_large- 1);
						if(found)
						out.emplace_back(large[pos]);
					}
				}
			};
	search(0, 0);
}

// without both LCA and early termination
template<indextype k = 3>
void intersect_hierarchical(keytype *small, indextype size_small,
		keytype *large, indextype size_large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	std::function<void(indextype depth, indextype rank)> search;
	search = [&]( indextype depth, indextype rank)-> void {

		// find the deepest common root
			for(indextype i = 0; i < k - 1; i++) {

				// now we can search the tree
				indextype pos = search_linearized_tree<k>(&found,
						small[PERFECT_SIZES[depth] + (k - 1) * rank + i],
						large, size_large - 1);
				if(found)
				out.emplace_back(large[pos]);

				// then the son nodes if exist
				if(PERFECT_SIZES[depth + 1] + (k * rank + i) * (k - 1) < size_small) {
					search(depth + 1,rank * k + i);
				}
			}

			if(PERFECT_SIZES[depth + 1] + (k * rank + k-1) * (k - 1) < size_small) {
				search(depth + 1,rank * k + k-1);
			}
		};
	search(0, 0);
}

/////////////////////////////////////////////////////////////////

/*** just a copy from sequential, need specialization before use ***/
template<indextype k>
void intersect_sequential_bounded(keytype *small, indextype size_small,
		keytype *large, indextype size_large, std::vector<keytype> & out) {
	out.clear();
	// linearly fetch
	int found = 0;

	for (indextype i = 0; i < size_small; i++) {
		/*indextype pos = */search_linearized_tree<k>(&found, small[i], large,
				size_large - 1);
//		std::cout << i << ": " << small[i] << ", " << large[pos] << std::endl;
		if (found)
			out.emplace_back(small[i]);
	}
}

template<>
void intersect_sequential_bounded<3>(keytype *small, indextype size_small,
		keytype *large, indextype size_large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	indextype height_large = getTreeHeight(size_large);
	indextype remainder_large = size_large - PERFECT_SIZES[height_large];

	indextype depth_rm, rank_rm;
	if (_LIKELY(remainder_large)) {
		depth_rm = height_large - 2;
		rank_rm = K_s[height_large - 2] - 1;
	} else {
		depth_rm = height_large - 1;
		rank_rm = K_s[height_large - 1] - 1;
	}

	indextype height_small = getTreeHeight(size_small);
	indextype remainder_small = size_small - PERFECT_SIZES[height_small];

	// here are two ranges used to store the boundries(depths and ranks)
	// for each level
	std::vector<indextype> range1((K_s[height_small - 1]) << 2);
	std::vector<indextype> range2((K_s[height_small - 1]) << 2);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = height_large - 1;
	p2[1] = p1[1] = 0;
	p1[2] = depth_rm;
	p1[3] = rank_rm - 1;

	// intersection begins here
	indextype depth = 0;
	if (_LIKELY(remainder_small))
		remainder_small = (size_small - PERFECT_SIZES[height_small - 1]) >> 1;
	else
		remainder_small = K_s[height_small - 1];

	for (indextype rank_limits = 1; depth < height_small - 1;
			depth++, rank_limits *= 3) {
		for (indextype rank = 0; rank < rank_limits; rank++) {
			/* search the first the first element */
			if (p1[rank << 1] == -1) {
				// note here we repeatedly set p2[rank*6] to -1
				// only aim to cover the position [0]
				// because [rank*6+6] cannot be [0]
				p2[rank * 6] = p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = p1[(rank << 1) + 2];
				p2[rank * 6 + 7] = p1[(rank << 1) + 3];
				continue;
			}
			indextype ii = 0;
			// -1 means no intersection exists here, we can safely skip
			while (p1[(rank << 1) + 2 + ii] == -1) {
				ii += 2;
			}

			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = p1[(rank << 1) + 2 + ii], _rrank = p1[(rank << 1)
							+ 3 + ii];

			_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: is it possible to use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[PERFECT_SIZES[depth] + (rank << 1)], large, _rdepth,
					_rrank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== PERFECT_SIZES[p1[rank << 1]]
									+ (p1[(rank << 1) + 1] << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6] = -1;
			} else if (_rdepth
					== p1[(rank << 1) + 2 + ii]&&
					_UNLIKELY(pos-PERFECT_SIZES[p1[(rank << 1) + 2 + ii]]-(p1[(rank << 1)
									+ 3 + ii]<<1)==2)) {
//			else if (_UNLIKELY(
//					pos - PERFECT_SIZES[_rdepth] - (_rrank << 1) == 2)) {

				// right end
				// its right sibling and right children are eliminated
				p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = p1[(rank << 1) + 2];
				p2[rank * 6 + 7] = p1[(rank << 1) + 3];
				continue;
			}
			_ldepth = p2[rank * 6 + 2] = _rdepth;
			_lrank = p2[rank * 6 + 3] = _rrank;
			_rdepth = p1[(rank << 1) + 2 + ii];
			_rrank = p1[(rank << 1) + 3 + ii];

			/* now the second element */
			_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[PERFECT_SIZES[depth] + (rank << 1) + 1], large,
					_rdepth, _rrank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== PERFECT_SIZES[p2[rank * 6 + 2]]
									+ (p2[rank * 6 + 3] << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6 + 2] = -1;
			} else if (_rdepth
					== p1[(rank << 1) + 2 + ii]&&
					_UNLIKELY(pos-PERFECT_SIZES[p1[(rank << 1) + 2 + ii]]-(p1[(rank << 1)
									+ 3 + ii]<<1)==2)) {
//			else if (_UNLIKELY(
//					pos - PERFECT_SIZES[_rdepth] - (_rrank << 1) == 2)) {
				// right end
				// its right children are eliminated
				p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = p1[(rank << 1) + 2];
				p2[rank * 6 + 7] = p1[(rank << 1) + 3];
				continue;
			}

			p2[rank * 6 + 4] = _rdepth;
			p2[rank * 6 + 5] = _rrank;
			p2[rank * 6 + 6] = p1[(rank << 1) + 2];
			p2[rank * 6 + 7] = p1[(rank << 1) + 3];
		}
		// swap p1 and p2
		std::swap(p1, p2);
	}

	/* handle the leaf level individually */
	for (indextype rank = 0; rank < remainder_small; rank++) {
		if (p1[rank << 1] == -1)
			continue;
		indextype ii = 0;
		while (p1[(rank << 1) + 2 + ii] == -1) {
			ii += 2;
		}

		/* search the first the first element */
		indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
				_rdepth = p1[(rank << 1) + 2 + ii], _rrank = p1[(rank << 1) + 3
						+ ii];

		_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: is it possible to use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[PERFECT_SIZES[depth] + (rank << 1)], large, _rdepth,
				_rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = p1[(rank << 1) + 2 + ii], _rrank = p1[(rank << 1) + 3 + ii];

		_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[PERFECT_SIZES[depth] + (rank << 1) + 1], large, _rdepth,
				_rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

// default k is 3 and no early termination included
void intersect_sequential_bounded_withoutET(keytype *small,
		indextype size_small, keytype *large, indextype size_large,
		std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	indextype height_large = getTreeHeight(size_large);
	indextype remainder_large = size_large - PERFECT_SIZES[height_large];

	indextype depth_rm, rank_rm;
	if (_LIKELY(remainder_large)) {
		depth_rm = height_large - 2;
		rank_rm = K_s[height_large - 2] - 1;
	} else {
		depth_rm = height_large - 1;
		rank_rm = K_s[height_large - 1] - 1;
	}

	indextype height_small = getTreeHeight(size_small);
	indextype remainder_small = size_small - PERFECT_SIZES[height_small];

	std::vector<indextype> range1(K_s[height_small - 1] << 2);
	std::vector<indextype> range2(K_s[height_small - 1] << 2);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

// leftmost node never changes
	p2[0] = p1[0] = height_large - 1;
	p2[1] = p1[1] = 0;
	p1[2] = depth_rm;
	p1[3] = rank_rm - 1;

// intersection begins here
	indextype depth = 0;
	if (_LIKELY(remainder_small))
		remainder_small = (size_small - PERFECT_SIZES[height_small - 1]) >> 1;
	else
		remainder_small = K_s[height_small - 1];

	for (indextype rank_limits = 1; depth < height_small - 1;
			depth++, rank_limits *= 3) {
		for (indextype rank = 0; rank < rank_limits; rank++) {
			/* search the first the first element */
			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

			_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: is it possible to use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[PERFECT_SIZES[depth] + (rank << 1)], large, _rdepth,
					_rrank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			_ldepth = p2[rank * 6 + 2] = _rdepth;
			_lrank = p2[rank * 6 + 3] = _rrank;
			_rdepth = p1[(rank << 1) + 2];
			_rrank = p1[(rank << 1) + 3];

			/* now the second element */
			_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[PERFECT_SIZES[depth] + (rank << 1) + 1], large,
					_rdepth, _rrank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			p2[rank * 6 + 4] = _rdepth;
			p2[rank * 6 + 5] = _rrank;
			p2[rank * 6 + 6] = p1[(rank << 1) + 2];
			p2[rank * 6 + 7] = p1[(rank << 1) + 3];
		}
		// swap p1 and p2
		std::swap(p1, p2);
	}

	/* handle the leaf level individually */
	for (indextype rank = 0; rank < remainder_small; rank++) {
		/* search the first the first element */
		indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
				_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

		_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: is it possible to use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[PERFECT_SIZES[depth] + (rank << 1)], large, _rdepth,
				_rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

		_lrank /= K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[PERFECT_SIZES[depth] + (rank << 1) + 1], large, _rdepth,
				_rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

/////////////////////////////////////////////////////////////////

// with LCA and early termination
template<indextype k = 3>
void intersect_hierarchical_bounded(keytype *small, indextype size_small,
		keytype *large, indextype size_large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;
// the rightmost node in the large tree
	indextype depth_rm, rank_rm;
	indextype height = getTreeHeight(size_large);
	indextype remainder = size_large - PERFECT_SIZES[height];
	if (_LIKELY(remainder)) {
		depth_rm = height - 2;
		rank_rm = K_s[height - 2] - 1;
	} else {
		depth_rm = height - 1;
		rank_rm = K_s[height - 1] - 1;
	}

	std::function<
			void(indextype depth, indextype rank, indextype ldepth,
					indextype lrank, indextype rdepth, indextype rrank)> search;
	search =
			[&]( indextype depth, indextype rank, indextype depth_left,
					indextype rank_left, indextype depth_right, indextype rank_right)-> void {

				indextype _lrank=rank_left,_ldepth=depth_left,_rrank,_rdepth;

				// find the deepest common root
				for(indextype i = 0; i < k - 1; i++) {
					_rrank=rank_right,_rdepth=depth_right;

//					if(_ldepth>_rdepth) {
//						_lrank /= K_s[_ldepth - _rdepth];
//					}
//					else if(_ldepth<_rdepth) {
//						_rrank /= K_s[_rdepth - _ldepth];
//					}
					_lrank/=K_s[_ldepth>_rdepth?_ldepth - _rdepth:0];
					_rrank/=K_s[_rdepth>_ldepth?_rdepth - _ldepth:0];
					_rdepth=min(_ldepth,_rdepth);

					while(_lrank != _rrank) {
						_rdepth--;
						_lrank /= k;
						_rrank /= k;
					}

					// now we can search the tree
					indextype pos = search_linearized_tree_bounded<k>(&found,
							small[PERFECT_SIZES[depth] + (k - 1) * rank + i],
							large, _rdepth, _rrank, size_large - 1);
					if(found)
					out.emplace_back(large[pos]);

//					std::cout<<"depth: "<<depth<<", rank: "<<(k - 1) * rank + i<<", value: "
//					<<small[PERFECT_SIZES[depth] + (k - 1) * rank + i]<<", found: "<<found<<std::endl;

					// then the left son nodes if exist
					// left early termination
					if(_LIKELY(!(pos==PERFECT_SIZES[depth_left]+rank_left*(k-1)) &&
									PERFECT_SIZES[depth + 1] + (k * rank + i) * (k - 1) < size_small)) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_rdepth,_rrank);
					}

					// do not modify again
					// right early termination
//					if(_UNLIKELY(pos-PERFECT_SIZES[_rdepth]-_rrank*(k-1)==k-1))
//					if(_UNLIKELY(pos-PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
//					if(_UNLIKELY(large[pos]>=large[PERFECT_SIZES[depth_right]+rank_right*(k-1)+k-2]))
					if(_rdepth==depth_right&&
							_UNLIKELY(pos-PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
					return;
					// update boundary nodes
					_ldepth=depth_left=_rdepth;
					_lrank=rank_left=_rrank;
				}

				if(PERFECT_SIZES[depth + 1] + (k * rank + k-1) * (k - 1) < size_small) {
					search(depth + 1,rank * k + k-1,depth_left,rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, height - 1, 0, depth_rm, rank_rm);
}

// with only LCA
template<indextype k = 3>
void intersect_hierarchical_bounded_withoutET(keytype *small,
		indextype size_small, keytype *large, indextype size_large,
		std::vector<keytype> & out) {
	out.clear();
	int found = 0;
// the rightmost node in the large tree
	indextype depth_rm, rank_rm;
	indextype height = getTreeHeight(size_large);
	indextype remainder = size_large - PERFECT_SIZES[height];
	if (_LIKELY(remainder)) {
		depth_rm = height - 2;
		rank_rm = K_s[height - 2] - 1;
	} else {
		depth_rm = height - 1;
		rank_rm = K_s[height - 1] - 1;
	}

	std::function<
			void(indextype depth, indextype rank, indextype ldepth,
					indextype lrank, indextype rdepth, indextype rrank)> search;
	search =
			[&]( indextype depth, indextype rank, indextype depth_left,
					indextype rank_left, indextype depth_right, indextype rank_right)-> void {

				indextype _lrank=rank_left,_ldepth=depth_left,_rrank,_rdepth;

				// find the deepest common root
				for(indextype i = 0; i < k - 1; i++) {
					_rrank=rank_right,_rdepth=depth_right;

//					if(_ldepth>_rdepth) {
//						_lrank /= K_s[_ldepth - _rdepth];
//					}
//					else if(_ldepth<_rdepth) {
//						_rrank /= K_s[_rdepth - _ldepth];
//					}
					_lrank/=K_s[_ldepth>_rdepth?_ldepth - _rdepth:0];
					_rrank/=K_s[_rdepth>_ldepth?_rdepth - _ldepth:0];
					_rdepth=min(_ldepth,_rdepth);

					while(_lrank != _rrank) {
						_rdepth--;
						_lrank /= k;
						_rrank /= k;
					}

					// now we can search the tree
					indextype pos = search_linearized_tree_bounded<k>(&found,
							small[PERFECT_SIZES[depth] + (k - 1) * rank + i],
							large, _rdepth, _rrank, size_large - 1);
					if(found)
					out.emplace_back(large[pos]);

//					std::cout<<"depth: "<<depth<<", rank: "<<(k - 1) * rank + i<<", value: "
//					<<small[PERFECT_SIZES[depth] + (k - 1) * rank + i]<<", found: "<<found<<std::endl;

					// then the left son nodes if exist
					// left early termination
					if(_LIKELY(PERFECT_SIZES[depth + 1] + (k * rank + i) * (k - 1) < size_small)) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_rdepth,_rrank);
					}

					// do not modify again
					// right early termination
//					if(_UNLIKELY(pos-PERFECT_SIZES[_rdepth]-_rrank*(k-1)==k-1))
//					if(_UNLIKELY(pos-PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
//					if(_UNLIKELY(large[pos]>=large[PERFECT_SIZES[depth_right]+rank_right*(k-1)+k-2]))
//					if(_rdepth==depth_right&&
//							_UNLIKELY(pos-PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
//					return;
					// update boundary nodes
					_ldepth=depth_left=_rdepth;
										_lrank=rank_left=_rrank;
				}

				if(PERFECT_SIZES[depth + 1] + (k * rank + k-1) * (k - 1) < size_small) {
					search(depth + 1,rank * k + k-1,depth_left,rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, height - 1, 0, depth_rm, rank_rm);
}

// with only early termination
// and its corresponding search
// algorithm will always start from root
template<indextype k = 3>
void intersect_hierarchical_bounded_withoutLCA(keytype *small,
		indextype size_small, keytype *large, indextype size_large,
		std::vector<keytype> & out) {
	out.clear();
	int found = 0;
// the rightmost node in the large tree
	indextype depth_rm, rank_rm;
	indextype height = getTreeHeight(size_large);
	indextype remainder = size_large - PERFECT_SIZES[height];
	if (_LIKELY(remainder)) {
		depth_rm = height - 2;
		rank_rm = K_s[height - 2] - 1;
	} else {
		depth_rm = height - 1;
		rank_rm = K_s[height - 1] - 1;
	}

	std::function<
			void(indextype depth, indextype rank, indextype ldepth,
					indextype lrank, indextype rdepth, indextype rrank)> search;
	search =
			[&]( indextype depth, indextype rank, indextype depth_left,
					indextype rank_left, indextype depth_right, indextype rank_right)-> void {

				// find the deepest common root
				for(indextype i = 0; i < k - 1; i++) {
					indextype _depth=0,_rank=0;

					// now we can search the tree
					indextype pos = search_linearized_tree_with_result<k>(&found,
							small[PERFECT_SIZES[depth] + (k - 1) * rank + i],
							large, _depth, _rank, size_large - 1);
					if(found)
					out.emplace_back(large[pos]);

					// then the son nodes if exist
					if(_LIKELY(!(pos==PERFECT_SIZES[depth_left]+rank_left*(k-1)) &&
									PERFECT_SIZES[depth + 1] + (k * rank + i) * (k - 1) < size_small)) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_depth,_rank);
					}
					if(_depth==depth_right&&
							_UNLIKELY(pos-PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
					return;
					// update boundary nodes
					depth_left = _depth;
					rank_left = _rank;
				}

				if(PERFECT_SIZES[depth + 1] + (k * rank + k-1) * (k - 1) < size_small) {
					search(depth + 1,rank * k + k-1,depth_left,rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, height - 1, 0, depth_rm, rank_rm);
}

}

