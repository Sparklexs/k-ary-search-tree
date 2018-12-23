/*
 * intersection_functions.hpp
 *
 *  Created on: 2018��12��16��
 *      Author: John
 */

#ifndef INCLUDE_INTERSECTION_FUNCTIONS_HPP_
#define INCLUDE_INTERSECTION_FUNCTIONS_HPP_

#include "k-ary_search_tree.hpp"
#include "search_functions.hpp"

namespace karytree {

///////////////////////////////////////////////////////////////////////////

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
// without both LCA and early termination
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

// with LCA and early termination
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
					// right early termination
					// if(_UNLIKELY(pos-PERFECT_SIZES[_rdepth]-_rrank*(k-1)==k-1))
					// if(_UNLIKELY(pos-PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
					// if(_UNLIKELY(large[pos]>=large[PERFECT_SIZES[depth_right]+rank_right*(k-1)+k-2]))
					// we may encounter one situation where the pos we found is exactly
					// the start position of next level and the given right boundary node
					// is the last node of the previous level,
					// then 'k-1' is not actually a sign of overstepping a node boundary
					// as it supposed to be. So it has to be ensured that we are at the
					// the same level with the right boundary nodes
					if(_rdepth==depth_right&&
							_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
					return;
					// update boundary nodes
					_ldepth=depth_left=_rdepth;
					_lrank=rank_left=_rrank;
				}

				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + k-1) * (k - 1) < small.size()) {
					search(depth + 1,rank * k + k-1,depth_left,
							rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, large.height() - 1, 0, large.rdepth(), large.rrank());
}

// with LCA and early termination
// get LCA via RMQ
template<indextype k>
void intersect_hierarchical_bounded_new(const ktree<k>& small,
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

					large.get_LCA(_ldepth,_lrank,_rdepth,_rrank);

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
					// right early termination
					// if(_UNLIKELY(pos-PERFECT_SIZES[_rdepth]-_rrank*(k-1)==k-1))
					// if(_UNLIKELY(pos-PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
					// if(_UNLIKELY(large[pos]>=large[PERFECT_SIZES[depth_right]+rank_right*(k-1)+k-2]))
					// the second sentence means we have
					if(_rdepth==depth_right&&
							_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[depth_right]-
									rank_right*(k-1)==k-1))
					return;
					// update boundary nodes
					_ldepth=depth_left=_rdepth;
					_lrank=rank_left=_rrank;
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

// invalid function (stupid copy from sequential)
// need specialization for reliable functionality
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
void intersect_sequential_bounded<3>(const ktree<3>& small,
		const ktree<3>& large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	std::vector<indextype> range1((ktree<3>::K_s[small.height() - 1]) << 1 + 1);
	std::vector<indextype> range2((ktree<3>::K_s[small.height() - 1]) << 1 + 1);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = large.height() - 1;
	p2[1] = p1[1] = 0;
	p1[2] = large.rdepth();
	p1[3] = large.rrank();

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
			indextype ii = 2;
			while (p1[(rank << 1) + ii] == -1) {
				ii += 2;
			}

			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = p1[(rank << 1) + ii], _rrank = p1[(rank << 1) + 1
							+ ii];
			indextype depth_right = _rdepth, rank_right = _rrank;

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
			} else if (depth_right
					== _rdepth&&_UNLIKELY(
							pos - ktree<3>::PERFECT_SIZES[depth_right] - (rank_right << 1)
							== 2)) {
				// right end
				// its right sibling and right children are eliminated
				p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}
			_ldepth = p2[rank * 6 + 2] = _rdepth;
			_lrank = p2[rank * 6 + 3] = _rrank;
			_rdepth = depth_right;
			_rrank = rank_right;

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
			} else if (depth_right == _rdepth&& _UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[depth_right]
					- (rank_right << 1)
					== 2)) {
				// right end
				// its right children are eliminated
				p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}

			p2[rank * 6 + 4] = _rdepth;
			p2[rank * 6 + 5] = _rrank;
			p2[rank * 6 + 6] = depth_right;
			p2[rank * 6 + 7] = rank_right;
		}
		// swap p1 and p2
		std::swap(p1, p2);
	}

	/* handle the leaf level individually */
	for (indextype rank = 0; rank <= small.leaf_rank(); rank++) {
		if (p1[rank << 1] == -1)
			continue;
		indextype ii = 2;
		while (p1[(rank << 1) + ii] == -1) {
			ii += 2;
		}

		/* search the first the first element */
		indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
				_rdepth = p1[(rank << 1) + ii], _rrank =
						p1[(rank << 1) + 1 + ii];
		indextype depth_right = _rdepth, rank_right = _rrank;

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

		// early right termination
		if (_rdepth == depth_right&&_UNLIKELY(
				pos - ktree<3>::PERFECT_SIZES[depth_right]
				- (rank_right << 1)
				== 2))
			continue;

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = depth_right, _rrank = rank_right;

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

void intersect_sequential_bounded_withoutLCA(const ktree<3>& small,
		const ktree<3>& large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	std::vector<indextype> range1((ktree<3>::K_s[small.height() - 1]) << 1 + 1);
	std::vector<indextype> range2((ktree<3>::K_s[small.height() - 1]) << 1 + 1);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = large.height() - 1;
	p2[1] = p1[1] = 0;
	p1[2] = large.rdepth();
	p1[3] = large.rrank();

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
			indextype ii = 2;
			while (p1[(rank << 1) + ii] == -1) {
				ii += 2;
			}

			indextype depth_left = p1[rank << 1], rank_left =
					p1[(rank << 1) + 1], depth_right = p1[(rank << 1) + ii],
					rank_right = p1[(rank << 1) + 1 + ii];
			indextype _depth = 0, _rank = 0;

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)],
					large.rootpointer(), _depth, _rank, large.size() - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== ktree<3>::PERFECT_SIZES[depth_left]
									+ (rank_left << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6] = -1;
			} else if (depth_right
					== _depth&&_UNLIKELY(
							pos - ktree<3>::PERFECT_SIZES[depth_right] - (rank_right << 1)
							== 2)) {
				// right end
				// its right sibling and right children are eliminated
				p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}
			depth_left = p2[rank * 6 + 2] = _depth;
			rank_left = p2[rank * 6 + 3] = _rank;
			_depth = 0;
			_rank = 0;

			/* now the second element */
			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
					large.rootpointer(), _depth, _rank, large.size() - 1);
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
			} else if (depth_right == _depth&& _UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[depth_right]
					- (rank_right << 1)
					== 2)) {
				// right end
				// its right children are eliminated
				p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}

			p2[rank * 6 + 4] = _depth;
			p2[rank * 6 + 5] = _rank;
			p2[rank * 6 + 6] = depth_right;
			p2[rank * 6 + 7] = rank_right;
		}
		// swap p1 and p2
		std::swap(p1, p2);
	}

	/* handle the leaf level individually */
	for (indextype rank = 0; rank <= small.leaf_rank(); rank++) {
		if (p1[rank << 1] == -1)
			continue;
		indextype ii = 2;
		while (p1[(rank << 1) + ii] == -1) {
			ii += 2;
		}

		/* search the first the first element */
		indextype depth_right = p1[(rank << 1) + ii], rank_right =
				p1[(rank << 1) + 1 + ii];
		indextype _depth = 0, _rank = 0;

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)],
				large.rootpointer(), _depth, _rank, large.size() - 1);
		if (found)
			out.emplace_back(large[pos]);

		// early right termination
		if (_depth == depth_right&&_UNLIKELY(
				pos - ktree<3>::PERFECT_SIZES[depth_right]
				- (rank_right << 1)
				== 2))
			continue;

		/* now the second element */
		_depth = 0, _rank = 0;

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
				large.rootpointer(), _depth, _rank, large.size() - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

// default k is 3 and no early termination is included
void intersect_sequential_bounded_withoutET(const ktree<3>& small,
		const ktree<3>& large, std::vector<keytype> & out) {
	out.clear();
	int found = 0;

	std::vector<indextype> range1(ktree<3>::K_s[small.height() - 1] << 1 + 1);
	std::vector<indextype> range2(ktree<3>::K_s[small.height() - 1] << 1 + 1);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = large.height() - 1;
	p2[1] = p1[1] = 0;
	p1[2] = large.rdepth();
	p1[3] = large.rrank();

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

// with only early termination
// and its corresponding search
// algorithm will always start from root
template<indextype k>
void intersect_hierarchical_bounded_withoutLCA(const ktree<k>& small,
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
					indextype pos = search_linearized_tree_bounded<k>(&found,
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
					if(_depth==depth_right&&
							_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[depth_right]-
									rank_right*(k-1)==k-1))
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

// with LCA and early termination
template<indextype k = 3>
void intersect_hierarchical_bounded_withoutET(const ktree<k>& small,
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

					// then the left son nodes if exist
					if(_LIKELY(ktree<k>::PERFECT_SIZES[depth + 1] +
									(k * rank + i) * (k - 1) < small.size())) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_rdepth,_rrank);
					}

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

}

namespace karytree_for_linearized_array {
using namespace karytree;

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
// without both LCA and early termination
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
			// the rightmost son node if exists
			if(ktree<k>::PERFECT_SIZES[depth + 1] +
					(k * rank + k-1) * (k - 1) < size_small) {
				search(depth + 1,rank * k + k-1);
			}
		};
	search(0, 0);
}

///////////////////////////////////////////////////////////////////////////

// with LCA and early termination
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
		rank_rm = ktree<k>::K_s[height - 2] - 1;
	} else {
		depth_rm = height - 1;
		rank_rm = ktree<k>::K_s[height - 1] - 1;
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
					if(_rdepth==depth_right&&
							_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[depth_right]-rank_right*(k-1)==k-1))
					return;
					// update boundary nodes
					_ldepth=depth_left=_rdepth;
					_lrank=rank_left=_rrank;
				}

				if(ktree<k>::PERFECT_SIZES[depth + 1] +
						(k * rank + k-1) * (k - 1) < size_small) {
					search(depth + 1,rank * k + k-1,depth_left,
							rank_left,depth_right,rank_right);
				}
			};
	search(0, 0, height - 1, 0, depth_rm, rank_rm);
}

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
		rank_rm = ktree<3>::K_s[height_large - 2] - 1;
	} else {
		depth_rm = height_large - 1;
		rank_rm = ktree<3>::K_s[height_large - 1] - 1;
	}

	indextype height_small = ktree<3>::getTreeHeight(size_small);
	indextype remainder_small = size_small
			- ktree<3>::PERFECT_SIZES[height_small];

	// for each node (not element) we prepare a pair of (depth, rank) of its LCA on large tree
	// and the leaf level has the most numerous nodes
	// actually we need n(nodes)+1 pairs
	std::vector<indextype> range1((ktree<3>::K_s[height_small - 1]) << 1 + 1);
	std::vector<indextype> range2((ktree<3>::K_s[height_small - 1]) << 1 + 1);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = height_large - 1;
	p2[1] = p1[1] = 0;
	p1[2] = depth_rm;
	p1[3] = rank_rm;

	// intersection begins here
	if (_LIKELY(remainder_small))
		remainder_small = (size_small
				- ktree<3>::PERFECT_SIZES[height_small - 1]) >> 1;
	else
		remainder_small = ktree<3>::K_s[height_small - 1];

	indextype depth = 0;
	for (indextype rank_limits = 1; depth < height_small - 1;
			depth++, rank_limits *= 3) {
		for (indextype rank = 0; rank < rank_limits; rank++) {
			/* search the first element */

			if (p1[rank << 1] == -1) {
				// note here we may repeatedly set p2[rank*6] to -1
				// only aim to cover the position [0]
				// because [rank*6+6] cannot be [0]
				p2[rank * 6] = p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = p1[(rank << 1) + 2];
				p2[rank * 6 + 7] = p1[(rank << 1) + 3];
				continue;
			}

			// skip over '-1'-points in the middle
			indextype ii = 2;
			while (p1[(rank << 1) + ii] == -1) {
				ii += 2;
			}
			indextype depth_right = p1[(rank << 1) + ii], rank_right = p1[(rank
					<< 1) + 1 + ii];
			indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
					_rdepth = depth_right, _rrank = rank_right;

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
			} else if (depth_right == _rdepth&& _UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[depth_right]
					- (rank_right << 1)
					== 2)) {
				// right end
				// its right sibling and right children are eliminated
				p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}
			// SXS:how to set [rank*6]?
			// we leave p2[rank*6] and p2[rank*6+1] unmodified
			// because they are set at the bottom of this for-loop
			_ldepth = p2[rank * 6 + 2] = _rdepth;
			_lrank = p2[rank * 6 + 3] = _rrank;
			_rdepth = depth_right;
			_rrank = rank_right;

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
			} else if (depth_right == _rdepth&&_UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[depth_right]
					- (rank_right << 1)
					== 2)) {
				// right end
				// its right children are eliminated
				p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}

			p2[rank * 6 + 4] = _rdepth;
			p2[rank * 6 + 5] = _rrank;
			p2[rank * 6 + 6] = depth_right;
			p2[rank * 6 + 7] = rank_right;
		}
		// swap p1 and p2
		std::swap(p1, p2);
	}

	/* handle the leaf level individually */
	for (indextype rank = 0; rank < remainder_small; rank++) {
		if (p1[rank << 1] == -1)
			continue;
		indextype ii = 2;
		while (p1[(rank << 1) + ii] == -1) {
			ii += 2;
		}

		/* search the first element */
		indextype _ldepth = p1[rank << 1], _lrank = p1[(rank << 1) + 1],
				_rdepth = p1[(rank << 1) + ii], _rrank =
						p1[(rank << 1) + 1 + ii];
		indextype depth_right = _rdepth, rank_right = _rrank;

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

		// early right termination
		if (_rdepth == depth_right&&_UNLIKELY(
				pos - ktree<3>::PERFECT_SIZES[depth_right]
				- (rank_right << 1)
				== 2))
			continue;

		/* now the second element */
		_ldepth = _rdepth;
		_lrank = _rrank;

		_rdepth = depth_right, _rrank = rank_right;

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

///////////////////////////////////////////////////////////////////////////

void intersect_sequential_bounded_withoutLCA(keytype *small,
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
		rank_rm = ktree<3>::K_s[height_large - 2] - 1;
	} else {
		depth_rm = height_large - 1;
		rank_rm = ktree<3>::K_s[height_large - 1] - 1;
	}

	indextype height_small = ktree<3>::getTreeHeight(size_small);
	indextype remainder_small = size_small
			- ktree<3>::PERFECT_SIZES[height_small];

	// for each node (not element) we prepare a pair of (depth, rank) of its LCA on large tree
	// and the leaf level has the most numerous nodes
	// actually we need n(nodes)+1 pairs
	std::vector<indextype> range1((ktree<3>::K_s[height_small - 1]) << 1 + 1);
	std::vector<indextype> range2((ktree<3>::K_s[height_small - 1]) << 1 + 1);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = height_large - 1;
	p2[1] = p1[1] = 0;
	p1[2] = depth_rm;
	p1[3] = rank_rm;

	// intersection begins here
	if (_LIKELY(remainder_small))
		remainder_small = (size_small
				- ktree<3>::PERFECT_SIZES[height_small - 1]) >> 1;
	else
		remainder_small = ktree<3>::K_s[height_small - 1];

	indextype depth = 0;
	for (indextype rank_limits = 1; depth < height_small - 1;
			depth++, rank_limits *= 3) {
		for (indextype rank = 0; rank < rank_limits; rank++) {
			/* search the first element */

			if (p1[rank << 1] == -1) {
				// note here we may repeatedly set p2[rank*6] to -1
				// only aim to cover the position [0]
				// because [rank*6+6] cannot be [0]
				p2[rank * 6] = p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = p1[(rank << 1) + 2];
				p2[rank * 6 + 7] = p1[(rank << 1) + 3];
				continue;
			}

			// skip over '-1'-points in the middle
			indextype ii = 2;
			while (p1[(rank << 1) + ii] == -1) {
				ii += 2;
			}
			indextype depth_right = p1[(rank << 1) + ii], rank_right = p1[(rank
					<< 1) + 1 + ii];
			indextype depth_left = p1[rank << 1], rank_left =
					p1[(rank << 1) + 1], _depth = 0, _rank = 0;

			// now we can search the tree
			indextype pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)], large,
					_depth, _rank, size_large - 1);
			if (found)
				out.emplace_back(large[pos]);

			// check early termination
			if (_UNLIKELY(
					pos
							== ktree<3>::PERFECT_SIZES[depth_left]
									+ (rank_left << 1))) {
				// left end
				// its left children are eliminated
				p2[rank * 6] = -1;
			} else if (depth_right == _depth&& _UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[depth_right]
					- (rank_right << 1)
					== 2)) {
				// right end
				// its right sibling and right children are eliminated
				p2[rank * 6 + 2] = p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}
			// SXS:how to set [rank*6]?
			// we leave p2[rank*6] and p2[rank*6+1] unmodified
			// because they are set at the bottom of this for-loop
			p2[rank * 6 + 2] = _depth;
			p2[rank * 6 + 3] = _rank;
			_depth = 0;
			_rank = 0;

			// now we can search the tree
			pos = search_linearized_tree_bounded<3>(&found,
					small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1],
					large, _depth, _rank, size_large - 1);
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
			} else if (depth_right == _depth&&_UNLIKELY(
					pos - ktree<3>::PERFECT_SIZES[depth_right]
					- (rank_right << 1)
					== 2)) {
				// right end
				// its right children are eliminated
				p2[rank * 6 + 4] = -1;
				p2[rank * 6 + 6] = depth_right;
				p2[rank * 6 + 7] = rank_right;
				continue;
			}

			p2[rank * 6 + 4] = _depth;
			p2[rank * 6 + 5] = _rank;
			p2[rank * 6 + 6] = depth_right;
			p2[rank * 6 + 7] = rank_right;
		}
		// swap p1 and p2
		std::swap(p1, p2);
	}

	/* handle the leaf level individually */
	for (indextype rank = 0; rank < remainder_small; rank++) {
		if (p1[rank << 1] == -1)
			continue;
		indextype ii = 2;
		while (p1[(rank << 1) + ii] == -1) {
			ii += 2;
		}

		/* search the first element */
		indextype depth_right = p1[(rank << 1) + ii], rank_right =
				p1[(rank << 1) + 1 + ii];
		indextype _depth = 0, _rank = 0;

		// now we can search the tree
		indextype pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1)], large,
				_depth, _rank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);

		// early right termination
		if (_depth == depth_right&&_UNLIKELY(
				pos - ktree<3>::PERFECT_SIZES[depth_right]
				- (rank_right << 1)
				== 2))
			continue;

		/* now the second element */
		_depth = 0, _rank = 0;

		// now we can search the tree
		pos = search_linearized_tree_bounded<3>(&found,
				small[ktree<3>::PERFECT_SIZES[depth] + (rank << 1) + 1], large,
				_depth, _rank, size_large - 1);
		if (found)
			out.emplace_back(large[pos]);
	}
}

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
		rank_rm = ktree<3>::K_s[height_large - 2] - 1;
	} else {
		depth_rm = height_large - 1;
		rank_rm = ktree<3>::K_s[height_large - 1] - 1;
	}

	indextype height_small = ktree<3>::getTreeHeight(size_small);
	indextype remainder_small = size_small
			- ktree<3>::PERFECT_SIZES[height_small];

	std::vector<indextype> range1(ktree<3>::K_s[height_small - 1] << 1 + 1);
	std::vector<indextype> range2(ktree<3>::K_s[height_small - 1] << 1 + 1);
	indextype *p1 = range1.data();
	indextype *p2 = range2.data();

	// leftmost node never changes
	p2[0] = p1[0] = height_large - 1;
	p2[1] = p1[1] = 0;
	p1[2] = depth_rm;
	p1[3] = rank_rm;

	// intersection begins here
	if (_LIKELY(remainder_small))
		remainder_small = (size_small
				- ktree<3>::PERFECT_SIZES[height_small - 1]) >> 1;
	else
		remainder_small = ktree<3>::K_s[height_small - 1];

	indextype depth = 0;
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

// with only early termination
// and its corresponding search
// algorithm will always start from root
template<indextype k>
void intersect_hierarchical_bounded_withoutLCA(keytype *small,
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
		rank_rm = ktree<k>::K_s[height - 2] - 1;
	} else {
		depth_rm = height - 1;
		rank_rm = ktree<k>::K_s[height - 1] - 1;
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
					indextype pos = search_linearized_tree_bounded<k>(&found,
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
					if(_depth==depth_right&&
							_UNLIKELY(pos-ktree<k>::PERFECT_SIZES[depth_right]-
									rank_right*(k-1)==k-1))
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

// with LCA and early termination
template<indextype k = 3>
void intersect_hierarchical_bounded_withoutET(keytype *small,
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
		rank_rm = ktree<k>::K_s[height - 2] - 1;
	} else {
		depth_rm = height - 1;
		rank_rm = ktree<k>::K_s[height - 1] - 1;
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

//					std::cout<<"depth: "<<depth<<", rank: "<<(k - 1) * rank + i<<", value: "
//					<<small[PERFECT_SIZES[depth] + (k - 1) * rank + i]<<", found: "<<found<<std::endl;

					// then the left son nodes if exist
					// left early termination
					if(_LIKELY(ktree<k>::PERFECT_SIZES[depth + 1] + (k * rank + i) * (k - 1) < size_small)) {
						search(depth + 1,rank * k + i,depth_left,rank_left,_rdepth,_rrank);
					}

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

}
#endif /* INCLUDE_INTERSECTION_FUNCTIONS_HPP_ */
