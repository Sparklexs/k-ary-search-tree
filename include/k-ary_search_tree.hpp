//============================================================================
// Name        : k-ary_search_tree.cpp
// Author      : Song
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================
#ifndef INCLUDE_K_ARY_SEARCH_TREE_HPP_
#define INCLUDE_K_ARY_SEARCH_TREE_HPP_

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

	// avoid signed integers
	while ((((vec[vec.size() - 1] + 1) * k - 1)
			& std::numeric_limits<indextype>::max())
			>> (sizeof(indextype) * 8 - 1) == ((indextype) 0)
			&& vec[vec.size() - 1] < (vec[vec.size() - 1] + 1) * k - 1) {
		vec.push_back((vec[vec.size() - 1] + 1) * k - 1);
	}

	return vec;
}

// k means the number of chunks in each level
template<indextype k>
inline std::vector<indextype> initializeKs() {
	std::vector<indextype> vec;
	vec.push_back(1);

	// avoid signed integers
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

	// number of elements in the leaf level
	// 0 if the tree is perfect
	indextype remainder() const {
		return m_remainder;
	}

	// number of nodes in the leaf level
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

	// number of nodes in the last level
	// last but one level if the tree is complete
	indextype rrank() const {
		return m_rrank;;
	}

	const keytype* rootpointer() const {
		return m_list.data();
	}

	void get_LCA(indextype leftdepth, indextype leftrank, indextype& rightdepth,
			indextype& rightrank) const {

		indextype l =
				m_occurrence[PERFECT_SIZES[leftdepth] / (k - 1) + leftrank];
		indextype r = m_occurrence[PERFECT_SIZES[rightdepth] / (k - 1)
				+ rightrank];

		if (l > r)
			std::swap(l, r);

		int h = 31 - __builtin_clz(static_cast<unsigned int>(r - l + 1));

		indextype rmq = get_min(m_st[l][h], m_st[r - (1 << h) + 1][h]);
		rightdepth = m_euler[rmq * 2];
		rightrank = m_euler[rmq * 2 + 1];

		/*********naive method to get LCA***********/
//			indextype ld = leftdepth, lr = leftrank, rd = rightdepth, rr =
//					rightrank;
//
//			leftrank /=
//					K_s[leftdepth > rightdepth ? leftdepth - rightdepth : 0];
//			rightrank /=
//					K_s[rightdepth > leftdepth ? rightdepth - leftdepth : 0];
//			rightdepth = std::min(leftdepth, rightdepth);
//
//			while (leftrank != rightrank) {
//				rightdepth--;
//				leftrank /= k;
//				rightrank /= k;
//			}
//		if (m_euler[rmq * 2] != rightdepth
//				|| m_euler[rmq * 2 + 1] != rightrank) {
//			std::cout << "(" << ld << "," << lr << ") (" << rd << "," << rr
//					<< "): " << std::flush;
//
//			std::cout << "(" << m_euler[rmq * 2] << "," << m_euler[rmq * 2 + 1]
//					<< ")    " << std::flush;
//
//			std::cout << "(" << rightdepth << "," << rightrank << ") "
//					<< std::endl;
//		}
	}

	/*
	 * @return height are always larger than 1
	 */
	static indextype getTreeHeight(indextype size) {
		indextype height = 0;
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
			m_leafrank = m_remainder / (k - 1) - 1;

			m_rdepth = m_height - 2;
			m_rrank = K_s[m_height < 2 ? 0 : m_height - 2] - 1;

		} else {
			m_rrank = m_leafrank = K_s[m_height - 1] - 1;
			m_rdepth = m_height - 1;
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

// TODO: there should be a series of functions to transform
// the index for a sorted array to its position in a tree and
// vice versa
// * index_to_pos
// * pos_to_index
// * next_index
// * next_pos
// however, these functions make no different to the search
// and intersection, they are presented only for functionally
// complete (maybe for range search)

private:
	std::vector<keytype> m_list;

	indextype m_height;

	indextype m_remainder;
	indextype m_leafrank;

	indextype m_ldepth, m_lrank;
	indextype m_rdepth, m_rrank;

// used for Lowest common ancestor
	std::vector<indextype> m_euler, m_depths, m_occurrence;
	std::vector<std::vector<indextype>> m_st;

	indextype get_min(indextype i, indextype j) const {
		return m_depths[i] < m_depths[j] ? i : j;
	}

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
		m_euler.push_back(depth);
		m_euler.push_back(rank);
		m_depths.push_back(depth);
		if (m_occurrence[PERFECT_SIZES[depth] / (k - 1) + rank]
				> m_depths.size() - 1)
			m_occurrence[PERFECT_SIZES[depth] / (k - 1) + rank] =
					m_depths.size() - 1;

		for (indextype i = 0; i < k - 1; i++) {
			index_tree[PERFECT_SIZES[depth] + (k - 1) * rank + i] =
					index_array[PERFECT_SIZES[height - 1] * (i + 1) + i];
		}

		if (--height) {
			// FIXME: for the leaf leaf level
			// we have only k-1 elements left
			// and they can be sequentially
			// copied to destination address
			for (indextype i = 0; i < k; i++) {
				reorder(index_array + i * (PERFECT_SIZES[height] + 1),
						rank * k + i, depth + 1, height, index_tree);
			}
		}

		m_euler.push_back(depth - 1);
		m_euler.push_back(rank / k);
		m_depths.push_back(depth - 1);
	}

// for now it only works for k equals 3
	void reorder_with_remainder(const keytype *index_array, indextype _size,
			indextype remainder, indextype rank, indextype depth,
			indextype height, keytype *index_tree) {
		m_euler.push_back(depth);
		m_euler.push_back(rank);
		m_depths.push_back(depth);

		if (m_occurrence[PERFECT_SIZES[depth] / (k - 1) + rank]
				> m_depths.size() - 1)
			m_occurrence[PERFECT_SIZES[depth] / (k - 1) + rank] =
					m_depths.size() - 1;

		indextype branch = remainder / (K_s[height - 2] * (k - 1));
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
			if (height != 2)
				reorder(
						right_sibling_start
								+ (i - branch) * PERFECT_SIZES[height - 2]
								+ (i + 1 - branch), rank * k + i + 1/*rank*/,
						depth + 1, height - 2, index_tree);
		}

		m_euler.push_back(depth - 1);
		m_euler.push_back(rank / k);
		m_depths.push_back(depth - 1);
	}

	void construct_tree_fast(const keytype* array, indextype _size) {
		m_st.assign(((_size >> 1) << 1) - 1, std::vector<indextype>());
		// large enough to be replaced
		m_occurrence.assign(_size / (k - 1), _size * k);

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

		// remove the invalid one
		m_depths.resize(m_depths.size() - 1);
		m_euler.resize(m_depths.size() * 2);

		// build sparse table for range minimum query
		for (indextype i = 0; i < m_depths.size(); i++)
			m_st[i].push_back(i);
		for (indextype j = 1; 1 << j <= m_depths.size(); j++)
			for (indextype i = 0; i + (1 << j) - 1 < m_depths.size(); i++)
				m_st[i].push_back(
						get_min(m_st[i][j - 1],
								m_st[i + (1 << (j - 1))][j - 1]));

//		for (indextype i = 0; i < m_depths.size(); i++)
//			m_st[i].push_back(get_min(i, std::min(i + 1, m_depths.size() - 1)));
//
//		for (indextype j = 1; (1 << j) <= m_depths.size(); j++)
//			for (indextype i = 0; i + (1 << j) + 1 < m_depths.size(); i++)
//				m_st[i].push_back(
//						get_min(m_st[i][j - 1],
//								m_st[std::min(m_depths.size() - 1, i + (1 << j))][j
//										- 1]));
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
}

#endif
