//============================================================================
// Name        : k-ary_search_tree.cpp
// Author      : Song
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <iostream>

#include "synthetic.h"
#include "timer.h"

namespace karytree {
typedef uint32_t keytype;
typedef size_t indextype;

template<indextype k>
inline std::vector<indextype> initializePerfectSizes() {
	std::vector<indextype> vec;
	vec.push_back(0);
	vec.push_back(k - 1);

	while ((((vec[vec.size() - 1] + 1) * k - 1)
			& std::numeric_limits<indextype>::max())
			>> (sizeof(indextype) * 8 - 1) == ((indextype) 0)
			&& vec[vec.size() - 1] < (vec[vec.size() - 1] + 1) * k - 1) {
		vec.push_back((vec[vec.size() - 1] + 1) * k - 1);
	}

	return vec;
}

template<indextype k>
inline std::vector<indextype> initializeKs() {
	std::vector<indextype> vec;
	vec.push_back(1);

	while ((vec[vec.size() - 1] * k & std::numeric_limits<indextype>::max())
			>> (sizeof(indextype) * 8 - 1) == ((indextype) 0)
			&& vec[vec.size() - 1] < vec[vec.size() - 1] * k) {
		vec.push_back(vec[vec.size() - 1] * k);
	}

	return vec;
}

template<indextype k>
class ktree {
public:
	indextype size() const {
		return m_list.size();
	}

	indextype height() const {
		return m_height;
	}

	indextype k_size() const {
		return k;
	}

	indextype remainder() const {
		return m_remainder;
	}

	indextype leaf_rank() const {
		return m_leafrank;
	}

	indextype ldepth() const {
		return m_ldepth;
	}

	indextype lrank() const {
		return m_lrank;
	}

	indextype rdepth() const {
		return m_rdepth;
	}

	indextype rrank() const {
		return m_rrank;;
	}

	const keytype* rootpointer() const {
		return m_list.data();
	}

	/*
	 * @return height are always larger than 1
	 */
	static indextype getTreeHeight(indextype size) {
		indextype height = 2;
		while (PERFECT_SIZES[height] < size)
			height++;
		return height;
	}

	ktree(const keytype* array, indextype _size) {
		m_list.resize(_size);
		m_height = getTreeHeight(_size);
		m_remainder = (_size) - PERFECT_SIZES[m_height];

		if (_LIKELY(m_remainder)) {
			m_remainder = (_size) - PERFECT_SIZES[m_height - 1];
			m_leafrank = m_remainder / (k - 1);

			m_rdepth = m_height - 2;
			m_rrank = K_s[m_height < 3 ? 0 : m_height - 3];

		} else {
			m_leafrank = K_s[m_height - 1];

			m_rdepth = m_height - 1;
			m_rrank = K_s[m_height - 2];

		}

		m_ldepth = m_height - 1;
		m_lrank = 0;

		construct_tree_fast(array, _size);
	}

	static const std::vector<indextype> PERFECT_SIZES;
	static const std::vector<indextype> K_s;
	inline const keytype& operator [](indextype pos) const {
		return m_list[pos];
	}
private:
	std::vector<keytype> m_list;

	indextype m_height;

	indextype m_remainder;
	indextype m_leafrank;

	indextype m_ldepth, m_lrank;
	indextype m_rdepth, m_rrank;

	void construct_tree(keytype* array, indextype _size) {
		indextype digits[100];
		m_list[PERFECT_SIZES[m_height - 1]] = array[0];

		// complete tree
		if (_LIKELY(m_remainder)) {
			indextype fringe_entry = (m_remainder - 1) * k / (k - 1);

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
				indextype depth = m_height - j - 1;
				indextype offset = 0;
				for (indextype p = pos - 1; p > j; p--) {
					offset = offset * k + digits[p] * (k - 1);
				}
				offset += digits[j] + PERFECT_SIZES[depth];
				m_list[offset] = array[i];
			}
			for (indextype i = fringe_entry + 1; i < _size; i++) {
				indextype pos = 0;
				indextype ind = i - m_remainder;
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
				indextype depth = m_height - j - 2;
				indextype offset = 0;
				for (indextype p = pos - 1 > pos ? 0 : pos - 1; p > j; p--) {
					offset = offset * k + digits[p] * (k - 1);
				}
				offset += digits[j] + PERFECT_SIZES[depth];
				m_list[offset] = array[i];
			}
		}
		// perfect tree
		else {
			for (indextype i = 1; i < _size; i++) {
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
				indextype depth = m_height - j - 1;
				indextype offset = 0;
				for (indextype p = pos - 1; p > j; p--) {
					offset = offset * k + digits[p] * (k - 1);
				}
				offset += digits[j] + PERFECT_SIZES[depth];
				m_list[offset] = array[i];
			}
		}
	}

	/*
	 * for nodes in the perfect tree
	 * @param index_array points to the start position of current subtree
	 * @param index_tree points to the very head of the root of the whole tree
	 */
	void reorder(const keytype *index_array, indextype rank, indextype depth,
			indextype height, keytype *index_tree) {
		if (height) {
			for (indextype i = 0; i < k - 1; i++) {
				index_tree[PERFECT_SIZES[depth] + (k - 1) * rank + i] =
						index_array[PERFECT_SIZES[height - 1] * (i + 1) + i];
			}
			if (--height) {
				for (indextype i = 0; i < k; i++) {
					reorder(index_array + i * (PERFECT_SIZES[height] + 1),
							rank * k + i, depth + 1, height, index_tree);
				}
			}
		}
	}

	// for now it only works for k equals 3
	void reorder_with_remainder(const keytype *index_array, indextype _size,
			indextype remainder, indextype rank, indextype depth,
			indextype height, keytype *index_tree) {
		indextype branch = remainder * k / (K_s[height - 1] * (k - 1));
		//if no complete tree exists, there will have a overlapping of one parent node

		// 1. the left sibling nodes are perfect trees with height of h-1
		for (indextype i = 0; i < branch; i++) {
			// here subtrees go before the parent nodes
			reorder(index_array + i * (PERFECT_SIZES[height - 1] + 1),
					rank * k + i/*rank*/, depth + 1/*depth*/, height - 1,
					index_tree);
			index_tree[PERFECT_SIZES[depth] + rank * (k - 1) + i] =
					index_array[PERFECT_SIZES[height - 1] * (i + 1) + i];
		}

		// 2. the 'branch'-th subtree is one complete tree with height of h-1
		// and it only exists when following condition is true
		if (remainder % ((K_s[height - 1] * (k - 1)) / k)) {
			indextype next_remainder = remainder
					% ((K_s[height - 1] * (k - 1)) / k);
			if (height == 2) {
				for (indextype i = 0; i < next_remainder; i++) {
					index_tree[_size - next_remainder + i] =
							index_array[PERFECT_SIZES[1] * branch + branch + i];
				}
			} else {
				reorder_with_remainder(
						index_array + branch * (PERFECT_SIZES[height - 1] + 1),
						_size, next_remainder, rank * k + branch, depth + 1,
						height - 1, index_tree);
			}
			branch++;
		}
		const keytype* right_sibling_start = index_array
				+ PERFECT_SIZES[height - 2] * branch + branch - 1 + remainder;
		// 3. the right sibling subtrees are perfect trees with height of h-2
		branch--;
		for (indextype i = branch; i < k - 1; i++) {
			// here is a little different from the first part, as the parent nodes
			// go before the subtrees, so rank*k + i '+ 1'
			index_tree[PERFECT_SIZES[depth] + rank * (k - 1) + i] =
					right_sibling_start[(i - branch) * PERFECT_SIZES[height - 2]
							+ (i - branch)];
			reorder(
					right_sibling_start
							+ (i - branch) * PERFECT_SIZES[height - 2]
							+ (i + 1 - branch), rank * k + i + 1/*rank*/,
					depth + 1, height - 2, index_tree);
		}
	}

	void construct_tree_fast(const keytype* array, indextype _size) {
		// complete tree
		if (_LIKELY(m_remainder)) {
			// remainder also means the number of elements in
			// the leaf level
			reorder_with_remainder(array, _size, m_remainder, 0, 0, m_height,
					m_list.data());
		}
		// perfect tree
		else {
			reorder(array, 0, 0, m_height, m_list.data());
		}
	}
};

template<>
const std::vector<indextype> ktree<3>::PERFECT_SIZES =
		initializePerfectSizes<3>();
template<>
const std::vector<indextype> ktree<5>::PERFECT_SIZES =
		initializePerfectSizes<5>();
template<>
const std::vector<indextype> ktree<9>::PERFECT_SIZES =
		initializePerfectSizes<9>();
template<>
const std::vector<indextype> ktree<17>::PERFECT_SIZES = initializePerfectSizes<
		17>();

template<>
const std::vector<indextype> ktree<3>::K_s = initializeKs<3>();
template<>
const std::vector<indextype> ktree<5>::K_s = initializeKs<5>();
template<>
const std::vector<indextype> ktree<9>::K_s = initializeKs<9>();
template<>
const std::vector<indextype> ktree<17>::K_s = initializeKs<17>();

///////////////////////////////////////////////////////////////////////////

/*
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
// always start search from root and return target (depth,rank)
template<indextype k>
indextype search_linearized_tree_bounded1(int *foundp, keytype goal,
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
		next = ktree<k>::PERFECT_SIZES[depth] + rank * (k - 1);
	} while (next < endindex);

	depth--;
	rank = (rank - m) / k;

	// m might be k -1, which means the largest value of the tree is less than goal
	*foundp = 0;
	return left + m;
}

///////////////////////////////////////////////////////////////////////////

template<indextype k>
void intersect_sequential(keytype *small, indextype size_small, keytype *large,
		indextype size_large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;
	for (indextype i = 0; i < size_small; i++) {
		search_linearized_tree<k>(&found, small[i], large, size_large - 1);
		if (found)
			out.emplace_back(small[i]);
	}
}

template<indextype k>
void intersect_sequential(const ktree<k>& small, const ktree<k>& large,
		std::vector<keytype> & out) {
	out.clear();

	int found = 0;
	for (indextype i = 0; i < small.size(); i++) {
		search_linearized_tree<k>(&found, small[i], large.rootpointer(),
				large.size() - 1);
		if (found)
			out.emplace_back(small[i]);
	}
}

///////////////////////////////////////////////////////////////////////////

// FIXME: assume each node is full of k-1 elements
template<indextype k>
void intersect_sorted(keytype *small, indextype size_small, keytype *large,
		indextype size_large, std::vector<keytype> & out) {
	out.clear();

	int found = 0;
	std::function<void(indextype rank, indextype depth)> search;
	search =
			[&]( indextype rank, indextype depth)-> void {
				for (indextype i = 0; i < k; i++) {
					if (ktree<k>::PERFECT_SIZES[depth + 1] +
							(k * rank + i) * (k - 1) < size_small) {
						search(rank * k + i, depth + 1);
					}
					if (i != k - 1) {
						indextype pos = search_linearized_tree<k>(&found,
								small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
								large, size_large- 1);
						if(found)
						out.emplace_back(large[pos]);
					}
				}
			};
	search(0, 0);
}

// FIXME: assume each node is full of k-1 elements
template<indextype k>
void intersect_sorted(const ktree<k>& small, const ktree<k>& large,
		std::vector<keytype> & out) {
	out.clear();

	int found = 0;
	std::function<void(indextype rank, indextype depth)> search;
	search =
			[&]( indextype rank, indextype depth)-> void {
				for (indextype i = 0; i < k; i++) {
					if (ktree<k>::PERFECT_SIZES[depth + 1] +
							(k * rank + i) * (k - 1) < small.size()) {
						search(rank * k + i, depth + 1);
					}
					if (i != k - 1) {
						indextype pos = search_linearized_tree<k>(&found,
								small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
								large.rootpointer(), large.size());
						if(found)
						out.emplace_back(large[pos]);
					}
				}
			};
	search(0, 0);
}

///////////////////////////////////////////////////////////////////////////

// FIXME: assume each node is full of k-1 elements
// without both DCR and early termination
template<indextype k>
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
						small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
						large, size_large - 1);
				if(found)
				out.emplace_back(large[pos]);

				// then the son nodes if exist
				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + i) * (k - 1) < size_small) {
					search(depth + 1,rank * k + i);
				}
			}

			if(ktree<k>::PERFECT_SIZES[depth + 1] +
					(k * rank + k-1) * (k - 1) < size_small) {
				search(depth + 1,rank * k + k-1);
			}
		};
	search(0, 0);
}

// FIXME: assume each node is full of k-1 elements
// without both DCR and early termination
template<indextype k>
void intersect_hierarchical(const ktree<k>& small, const ktree<k>& large,
		std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	std::function<void(indextype depth, indextype rank)> search;
	search = [&]( indextype depth, indextype rank)-> void {

		// find the deepest common root
			for(indextype i = 0; i < k - 1; i++) {

				// now we can search the tree
				indextype pos = search_linearized_tree<k>(&found,
						small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
						large.rootpointer(), large.size());
				if(found)
				out.emplace_back(large[pos]);

				// then the son nodes if exist
				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + i) * (k - 1) < small.size()) {
					search(depth + 1,rank * k + i);
				}
			}

			if(ktree<k>::PERFECT_SIZES[depth + 1] +
					(k * rank + k-1) * (k - 1) < small.size()) {
				search(depth + 1,rank * k + k-1);
			}
		};
	search(0, 0);
}

///////////////////////////////////////////////////////////////////////////

// invalid function (stupid copy from sequential)
// need specialization for reliable functionality
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

template<indextype k>
void intersect_sequential_bounded(const ktree<k>& small, const ktree<k>& large,
		std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	for (indextype i = 0; i < small.size(); i++) {
		search_linearized_tree<k>(&found, small[i], large.rootpointer(),
				large.size() - 1);
		if (found)
			out.emplace_back(small[i]);
	}
}

///////////////////////////////////////////////////////////////////////////

template<>
void intersect_sequential_bounded<3>(keytype *small, indextype size_small,
		keytype *large, indextype size_large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	indextype height_large = ktree<3>::getTreeHeight(size_large);
	indextype remainder_large = size_large
			- ktree<3>::PERFECT_SIZES[height_large];

	indextype depth_rm, rank_rm;
	if (_LIKELY(remainder_large)) {
		depth_rm = height_large - 2;
		rank_rm = ktree<3>::K_s[height_large - 2];
	} else {
		depth_rm = height_large - 1;
		rank_rm = ktree<3>::K_s[height_large - 1];
	}

	indextype height_small = ktree<3>::getTreeHeight(size_small);
	indextype remainder_small = size_small
			- ktree<3>::PERFECT_SIZES[height_small];

	std::vector<indextype> range1((ktree<3>::K_s[height_small - 1]) << 2);
	std::vector<indextype> range2((ktree<3>::K_s[height_small - 1]) << 2);
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
		remainder_small = (size_small
				- ktree<3>::PERFECT_SIZES[height_small - 1]) >> 1;
	else
		remainder_small = ktree<3>::K_s[height_small - 1];

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
			while (p1[(rank << 1) + 2 + ii] == -1) {
				ii += 2;
			}

			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = p1[(rank << 1) + 2 + ii], _rrank = p1[(rank << 1)
							+ 3 + ii];

			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: is it possible to use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)], large,
					_rdepth, _rrank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== ktree<3>::PERFECT_SIZES[p1[rank << 1]]
									+ (p1[(rank << 1) + 1] << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6] = -1;
			} else if (_UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[_rdepth] - (_rrank << 1)
							== 2)) {
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
			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
					large, _rdepth, _rrank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== ktree<3>::PERFECT_SIZES[p2[rank * 6 + 2]]
									+ (p2[rank * 6 + 3] << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6 + 2] = -1;
			} else if (_UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[_rdepth] - (_rrank << 1)
							== 2)) {
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

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: is it possible to use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)], large,
				_rdepth, _rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = p1[(rank << 1) + 2 + ii], _rrank = p1[(rank << 1) + 3 + ii];

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1], large,
				_rdepth, _rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

template<>
void intersect_sequential_bounded<3>(const ktree<3>& small,
		const ktree<3>& large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	std::vector<indextype> range1((ktree<3>::K_s[small.height() - 1]) << 2);
	std::vector<indextype> range2((ktree<3>::K_s[small.height() - 1]) << 2);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = large.height() - 1;
	p2[1] = p1[1] = 0;
	p1[2] = large.rdepth();
	p1[3] = large.rrank() - 1;

	// intersection begins here
	indextype depth = 0;

	for (indextype rank_limits = 1; depth < small.height() - 1;
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
			while (p1[(rank << 1) + 2 + ii] == -1) {
				ii += 2;
			}

			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = p1[(rank << 1) + 2 + ii], _rrank = p1[(rank << 1)
							+ 3 + ii];

			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: is it possible to use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)],
					large.rootpointer(), _rdepth, _rrank, large.size() - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== ktree<3>::PERFECT_SIZES[p1[rank << 1]]
									+ (p1[(rank << 1) + 1] << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6] = -1;
			} else if (_UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[_rdepth] - (_rrank << 1)
							== 2)) {
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
			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
					large.rootpointer(), _rdepth, _rrank, large.size() - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== ktree<3>::PERFECT_SIZES[p2[rank * 6 + 2]]
									+ (p2[rank * 6 + 3] << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6 + 2] = -1;
			} else if (_UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[_rdepth] - (_rrank << 1)
							== 2)) {
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
	for (indextype rank = 0; rank < small.remainder(); rank++) {
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

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: is it possible to use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)],
				large.rootpointer(), _rdepth, _rrank, large.size() - 1);
		if (found)
			out.emplace_back(large[pos]);

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = p1[(rank << 1) + 2 + ii], _rrank = p1[(rank << 1) + 3 + ii];

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
				large.rootpointer(), _rdepth, _rrank, large.size() - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

///////////////////////////////////////////////////////////////////////////

// default k is 3 and no early termination is included
void intersect_sequential_bounded_withoutET(keytype *small,
		indextype size_small, keytype *large, indextype size_large,
		std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	indextype height_large = ktree<3>::getTreeHeight(size_large);
	indextype remainder_large = size_large
			- ktree<3>::PERFECT_SIZES[height_large];

	indextype depth_rm, rank_rm;
	if (_LIKELY(remainder_large)) {
		depth_rm = height_large - 2;
		rank_rm = ktree<3>::K_s[height_large - 2];
	} else {
		depth_rm = height_large - 1;
		rank_rm = ktree<3>::K_s[height_large - 1];
	}

	indextype height_small = ktree<3>::getTreeHeight(size_small);
	indextype remainder_small = size_small
			- ktree<3>::PERFECT_SIZES[height_small];

	std::vector<indextype> range1(ktree<3>::K_s[height_small - 1] << 2);
	std::vector<indextype> range2(ktree<3>::K_s[height_small - 1] << 2);
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
		remainder_small = (size_small
				- ktree<3>::PERFECT_SIZES[height_small - 1]) >> 1;
	else
		remainder_small = ktree<3>::K_s[height_small - 1];

	for (indextype rank_limits = 1; depth < height_small - 1;
			depth++, rank_limits *= 3) {
		for (indextype rank = 0; rank < rank_limits; rank++) {
			/* search the first the first element */
			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: is it possible to use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)], large,
					_rdepth, _rrank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			_ldepth = p2[rank * 6 + 2] = _rdepth;
			_lrank = p2[rank * 6 + 3] = _rrank;
			_rdepth = p1[(rank << 1) + 2];
			_rrank = p1[(rank << 1) + 3];

			/* now the second element */
			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
					large, _rdepth, _rrank, size_large - 1);
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

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: is it possible to use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)], large,
				_rdepth, _rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1], large,
				_rdepth, _rrank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

// default k is 3 and no early termination is included
void intersect_sequential_bounded_withoutET(const ktree<3>& small,
		const ktree<3>& large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	std::vector<indextype> range1(ktree<3>::K_s[small.height() - 1] << 2);
	std::vector<indextype> range2(ktree<3>::K_s[small.height() - 1] << 2);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = large.height() - 1;
	p2[1] = p1[1] = 0;
	p1[2] = large.rdepth();
	p1[3] = large.rrank() - 1;

	// intersection begins here
	indextype depth = 0;

	for (indextype rank_limits = 1; depth < small.height() - 1;
			depth++, rank_limits *= 3) {
		for (indextype rank = 0; rank < rank_limits; rank++) {
			/* search the first the first element */
			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: is it possible to use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)],
					large.rootpointer(), _rdepth, _rrank, large.size() - 1);
			if (found)
				out.emplace_back(large[pos]);

			_ldepth = p2[rank * 6 + 2] = _rdepth;
			_lrank = p2[rank * 6 + 3] = _rrank;
			_rdepth = p1[(rank << 1) + 2];
			_rrank = p1[(rank << 1) + 3];

			/* now the second element */
			_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
			_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
			_rdepth = min(_ldepth, _rdepth);

			// TODO: use xor to replace this while-loop
			while (_lrank != _rrank) {
				_rdepth--;
				_lrank /= 3;
				_rrank /= 3;
			}

			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
					large.rootpointer(), _rdepth, _rrank, large.size() - 1);
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
	for (indextype rank = 0; rank < small.remainder(); rank++) {
		/* search the first the first element */
		indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
				_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: is it possible to use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)],
				large.rootpointer(), _rdepth, _rrank, large.size() - 1);
		if (found)
			out.emplace_back(large[pos]);

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = p1[(rank << 1) + 2], _rrank = p1[(rank << 1) + 3];

		_lrank /= ktree<3>::K_s[_ldepth > _rdepth ? _ldepth - _rdepth : 0];
		_rrank /= ktree<3>::K_s[_rdepth > _ldepth ? _rdepth - _ldepth : 0];
		_rdepth = min(_ldepth, _rdepth);

		// TODO: use xor to replace this while-loop
		while (_lrank != _rrank) {
			_rdepth--;
			_lrank /= 3;
			_rrank /= 3;
		}

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
				large.rootpointer(), _rdepth, _rrank, large.size() - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

///////////////////////////////////////////////////////////////////////////

// with DCR and early termination
template<indextype k>
void intersect_hierarchical_bounded(keytype *small, indextype size_small,
		keytype *large, indextype size_large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;
	// the rightmost node in the large tree
	indextype depth_rm, rank_rm;
	indextype height = ktree<k>::getTreeHeight(size_large);
	indextype remainder = size_large - ktree<k>::PERFECT_SIZES[height];
	if (_LIKELY(remainder)) {
		depth_rm = height - 2;
		rank_rm = ktree<k>::K_s[height - 3];
	} else {
		depth_rm = height - 1;
		rank_rm = ktree<k>::K_s[height - 2];
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
					_lrank/=ktree<k>::K_s[_ldepth>_rdepth?_ldepth - _rdepth:0];
					_rrank/=ktree<k>::K_s[_rdepth>_ldepth?_rdepth - _ldepth:0];
					_rdepth=min(_ldepth,_rdepth);

					while(_lrank != _rrank) {
						_rdepth--;
						_lrank /= k;
						_rrank /= k;
					}

					// now we can search the tree
					indextype pos = search_linearized_tree_bounded<k>(&found,
							small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
							large, _rdepth, _rrank, size_large - 1);
					if(found)
					out.emplace_back(large[pos]);

					// then the son nodes if exist
					if(_LIKELY(!(pos==ktree<k>::PERFECT_SIZES[depth_left]+rank_left*(k-1)) &&
									ktree<k>::PERFECT_SIZES[depth + 1] +
									(k * rank + i) * (k - 1) < size_small)) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_rdepth,_rrank);
					}

					// do not modify again
					if(_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[_rdepth]-_rrank*(k-1)==k-1))
					return;
					// update boundary nodes
					depth_left=_rdepth;
					rank_left=_rrank;
				}

				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + k-1) * (k - 1) < size_small) {
					search(depth + 1,rank * k + k-1,depth_left,
							rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, height - 1, 0, depth_rm, rank_rm);
}

// with DCR and early termination
template<indextype k>
void intersect_hierarchical_bounded(const ktree<k>& small,
		const ktree<k>& large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

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

					_lrank/=ktree<k>::K_s[_ldepth>_rdepth?_ldepth - _rdepth:0];
					_rrank/=ktree<k>::K_s[_rdepth>_ldepth?_rdepth - _ldepth:0];
					_rdepth=min(_ldepth,_rdepth);

					while(_lrank != _rrank) {
						_rdepth--;
						_lrank /= k;
						_rrank /= k;
					}

					// now we can search the tree
					indextype pos = search_linearized_tree_bounded<k>(&found,
							small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
							large.rootpointer(), _rdepth, _rrank, large.size() - 1);
					if(found)
					out.emplace_back(large[pos]);

					// then the son nodes if exist
					if(_LIKELY(!(pos==ktree<k>::PERFECT_SIZES[depth_left]+rank_left*(k-1)) &&
									ktree<k>::PERFECT_SIZES[depth + 1] +
									(k * rank + i) * (k - 1) < small.size())) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_rdepth,_rrank);
					}

					// do not modify again
					if(_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[_rdepth]-_rrank*(k-1)==k-1))
					return;
					// update boundary nodes
					depth_left=_rdepth;
					rank_left=_rrank;
				}

				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + k-1) * (k - 1) < small.size()) {
					search(depth + 1,rank * k + k-1,depth_left,
							rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, large.height() - 1, 0, large.rdepth(), large.rrank());
}

///////////////////////////////////////////////////////////////////////////

// with only early termination
// and its corresponding search
// algorithm will always start from root
template<indextype k>
void intersect_hierarchical_bounded_withoutDCR(keytype *small,
		indextype size_small, keytype *large, indextype size_large,
		std::vector<keytype> & out) {
	out.clear();
	int found = 0;
	// the rightmost node in the large tree
	indextype depth_rm, rank_rm;
	indextype height = ktree<k>::getTreeHeight(size_large);
	indextype remainder = size_large - ktree<k>::PERFECT_SIZES[height];
	if (_LIKELY(remainder)) {
		depth_rm = height - 2;
		rank_rm = ktree<k>::K_s[height - 3];
	} else {
		depth_rm = height - 1;
		rank_rm = ktree<k>::K_s[height - 2];
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
					indextype pos = search_linearized_tree_bounded1<k>(&found,
							small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
							large, _depth, _rank, size_large - 1);
					if(found)
					out.emplace_back(large[pos]);

					// then the son nodes if exist
					if(_LIKELY(!(pos==ktree<k>::PERFECT_SIZES[depth_left]+rank_left*(k-1)) &&
									ktree<k>::PERFECT_SIZES[depth + 1] +
									(k * rank + i) * (k - 1) < size_small)) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_depth,_rank);
					}
					if(_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[_depth]-_rank*(k-1)==k-1))
					return;
					// update boundary nodes
					depth_left = _depth;
					rank_left = _rank;
				}

				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + k-1) * (k - 1) < size_small) {
					search(depth + 1,rank * k + k-1,depth_left,
							rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, height - 1, 0, depth_rm, rank_rm);
}

// with only early termination
// and its corresponding search
// algorithm will always start from root
template<indextype k>
void intersect_hierarchical_bounded_withoutDCR(const ktree<k>& small,
		const ktree<k>& large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

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
					indextype pos = search_linearized_tree_bounded1<k>(&found,
							small[ktree<k>::PERFECT_SIZES[depth] + (k - 1) * rank + i],
							large.rootpointer(), _depth, _rank, large.size() - 1);
					if(found)
					out.emplace_back(large[pos]);

					// then the son nodes if exist
					if(_LIKELY(!(pos==ktree<k>::PERFECT_SIZES[depth_left]+rank_left*(k-1)) &&
									ktree<k>::PERFECT_SIZES[depth + 1] +
									(k * rank + i) * (k - 1) < small.size())) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_depth,_rank);
					}
					if(_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[_depth]-_rank*(k-1)==k-1))
					return;
					// update boundary nodes
					depth_left = _depth;
					rank_left = _rank;
				}

				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + k-1) * (k - 1) < small.size()) {
					search(depth + 1,rank * k + k-1,depth_left,
							rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, large.height() - 1, 0, large.rdepth(), large.rrank());
}
}
