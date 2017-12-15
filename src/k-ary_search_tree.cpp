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

namespace karytree {
typedef uint32_t keytype;
typedef size_t indextype;

const static indextype PERFECT_SIZES[] = { 0, 2, 8, 26, 80, 242, 728, 2186,
		6560, 19682, 59048, 177146, 531440, 1594322, 4782968, 14348906,
		43046720, 129140162, 387420488, 1162261466, 3486784400 /*, 10460353202, 31381059608,
 94143178826, 282429536480, 847288609442, 2541865828328, 7625597484986,
 22876792454960, 68630377364882, 205891132094648, 617673396283946,
 1853020188851840, 5559060566555522*/};

// the static array PERFECT_SIZES is only used for k equals 3
// so we have to calculate actual array lengths when k varies.
template<uint32_t k = 3>
inline void initializeArrayLength(std::vector<keytype>& vec) {
	vec[0] = k - 1;
	for (indextype i = 1; i < vec.size(); i++) {
		vec[i] = (vec[i - 1] + 1) * k - 1;
	}
}

/*
 * @return height are always larger than 1
 */
indextype getTreeHeight(indextype size) {
	indextype height = 2;
	while (PERFECT_SIZES[height] < size)
		height++;
	return height;
}

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

// this function should be moved into the ktree class
template<indextype k = 3>
indextype search_linearized_tree(keytype goal, const keytype *target,
		indextype endindex) {
	indextype left = 0, next = 0;
	indextype m;
	while (next < endindex) {
		left = next;
		for (m = 0; m < k - 1; m++) {
			if (target[left + m] == goal)
				return left + m;
			if (target[left + m] > goal)
				break;
		}
		next = (left + 1) * k + m * (k - 1) - 1;
	}
	return left + m;
}

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
			offset_next = offset / (k - 1)
					* (PERFECT_SIZES[height - depth]
							- PERFECT_SIZES[height - depth - 1])
					+ (offset % (k - 1) + 1)
							* (PERFECT_SIZES[height - depth - 1]
									- PERFECT_SIZES[height - depth - 2]);
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
			offset_next = offset / (k - 1)
					* (PERFECT_SIZES[height - depth]
							- PERFECT_SIZES[height - depth - 1])
					+ (offset % (k - 1) + 1)
							* (PERFECT_SIZES[height - depth - 1]
									- PERFECT_SIZES[height - depth - 2]);
			return offset_next + PERFECT_SIZES[depth_next];
		}
		// only elements from leaf nodes, their next elements can be
		// right next to them or grow into upper nodes
		LEAF: if ((offset + 1) % (k - 1) != 0) {
			return pos + 1;
		} else if (offset
				== PERFECT_SIZES[height] - PERFECT_SIZES[height - 1] - 1) {
			// manually circular-join the last leaf to the root
			return 0;
		} else {
			ASCEND: indextype sum_of_middle_points = 0;
			while ((offset + 1)
					% (PERFECT_SIZES[height - depth_next - 1]
							- PERFECT_SIZES[height - depth_next - 2]) != 0) {
				sum_of_middle_points += (offset + 1)
						/ (PERFECT_SIZES[height - depth_next - 1]
								- PERFECT_SIZES[height - depth_next - 2]);
				depth_next++;
			}
			offset_next = (offset + 1)
					/ (PERFECT_SIZES[height - depth_next - 1]
							- PERFECT_SIZES[height - depth_next - 2])
					- sum_of_middle_points - 1;
			return offset_next + PERFECT_SIZES[depth_next];
		}
	}
}

template<indextype k = 3>
void construct_tree(keytype* array, indextype size) {
	indextype height = getTreeHeight(size);
	indextype remainder = (size) - PERFECT_SIZES[height];
	keytype tmp[size];
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
			while (digits[j] == k - 1)		// == 2?
				j++;
			indextype depth = height - j - 2;
			indextype offset = 0;
			for (indextype p = pos - 1; p > j; p--) {
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
			while (digits[j] == k - 1)		// == 2?
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
	memmove(array, tmp, size * sizeof(keytype));
}

template<indextype k = 3>
class ktree {
public:
	indextype size() {
		return m_size;
	}

	indextype height() {
		return m_height;
	}

	indextype k_size() {
		return k;
	}

	ktree(const keytype* array, indextype _size) {
		m_size = _size;
		m_list.resize(m_size);

		m_height = getTreeHeight(_size);
		indextype remainder = (_size) - PERFECT_SIZES[m_height];

		// complete tree
		if (_LIKELY(remainder)) {
			remainder = (_size) - PERFECT_SIZES[m_height - 1];
			indextype fringe_entry = (remainder - 1) * k / (k - 1);

			indextype digits[100];
			m_list[PERFECT_SIZES[m_height - 1]] = array[0];
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
				indextype depth = m_height - j - 2;
				indextype offset = 0;
				for (indextype p = pos - 1; p > j; p--) {
					offset = offset * k + digits[p] * (k - 1);
				}
				offset += digits[j] + PERFECT_SIZES[depth];
				m_list[offset] = array[i];
			}
		}
		// perfect tree
		else {
			indextype digits[100];
			m_list[PERFECT_SIZES[m_height - 1]] = array[0];
			for (indextype i = 1; i < _size; i++) {
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
		}
	}

	class sorted_iterator;
	friend class sorted_iterator;
	sorted_iterator so_begin() const {
		return sorted_iterator(this, PERFECT_SIZES[m_height - 1]);
	}
	sorted_iterator so_end() const {
		return sorted_iterator(this, m_size - 1);
	}
	std::vector<keytype>::iterator se_begin() const {
		return m_list.begin();
	}
	std::vector<keytype>::iterator se_end() const {
		return m_list.end();
	}
	class sorted_iterator: public std::iterator<std::forward_iterator_tag,
			keytype> {
	public:
		keytype const& operator*() const {
			return m_tree->m_list[m_pos];
		}
		keytype const& operator->() const {
			return m_tree->m_list[m_pos];
		}
		iterator& operator++() {

			return *this;
		}
		bool operator==(sorted_iterator const& other) const {
			assert(m_tree == other.m_tree);
			return m_pos == other.m_pos;
		}
		bool operator!=(sorted_iterator const& other) const {
			return !(*this == other);
		}

	private:
		friend class ktree;
		sorted_iterator(ktree const* _tree, indextype pos) :
				m_tree(_tree), m_pos(pos) {

		}
		ktree* m_tree;
		indextype m_pos;
	};
private:
	indextype m_size;
	indextype m_height;
	std::vector<keytype> m_list;
};

}

int main(void) {
	using namespace karytree;

	std::vector<keytype> vec(32);
	initializeArrayLength(vec);

	ClusteredDataGenerator sdg;
	std::vector<keytype> ans = sdg.generate(20, 20);
	construct_tree(ans.data(), ans.size());

	for (indextype i = 0; i < 20; i++) {
		indextype ind = next_pos(i, 20);
		std::cout << i << ": " << ind << std::endl;
	}

	keytype goal = 1;
	for (uint32_t goal = 0; goal < 26; goal++) {
		indextype pos = search_sorted_array(goal, ans.data(), ans.size() - 1);
		if (ans[pos] != goal)
			std::cout << "bad! Found ";
		else
			std::cout << "good! Found ";
		std::cout << ans[pos] << " for " << goal << std::endl;
	}
	return 0;
}
