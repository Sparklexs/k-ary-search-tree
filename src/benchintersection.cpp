/*
 * benchintersection.cpp
 *
 *  Created on: 2018��1��15��
 *      Author: John
 */

#include "k-ary_search_tree_trial.hpp"
#include "intersection_functions.hpp"

template<karytree::indextype k = 3>
void run() {
	using namespace karytree;

	vector<uint32_t> out;

	uint32_t MaxBit = 31; // largest bit-length of element
	uint32_t minlength;
	size_t REPETITION = 1;
	size_t CASES = 20;

	WallClockTimer timer;
	vector<float> intersectionsratios = { 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
			0.70, 0.80, 0.90, 1.00 };
	vector<uint32_t> sizeratios = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
			2048, 4096, 8192 };

	std::map<std::string, size_t> times;
	for (uint32_t msb = 10; msb <= 14; ++msb) {
		minlength = 1U << msb;

		for (float ir : intersectionsratios) {
			printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
			for (uint32_t sr : sizeratios) {
				printf("  size ratio: \e[32m%4d\e[0m\n", sr);

				if (sr > 1000)
					REPETITION = 100;
				else
					REPETITION = 200;

				for (uint32_t num = 2; num < 3; num++) {
					times["hierarchical"] = 0;
					times["sequential"] = 0;
					times["sorted"] = 0;

					times["hierarchical_bounded"] = 0;
					times["hierarchical_bounded_noLCA"] = 0;
					times["sequential_bounded"] = 0;
					times["sequential_bounded_noET"] = 0;
					time_t t = time(nullptr);
					tm* _tm = localtime(&t);
					printf("%02d:%02d:%02d> num: \e[32m%2d\e[0m  \n",
							_tm->tm_hour, _tm->tm_min, _tm->tm_sec, num);

					for (uint32_t instance = 0; instance < CASES; instance++) {
						mySet multiset;

						ClusteredDataGenerator cdg;
//						std::cout << "we are here!" << std::endl;
//						return;
						multiset = genMultipleSets(cdg, minlength, num,
								1U << MaxBit, static_cast<float>(sr), ir);

						// build k-ary tree
						auto it = multiset.begin();
						std::vector<keytype> small(std::move(*it++));
						std::vector<keytype> large(std::move(*it));

						ktree<k> smalltree(small.data(), small.size());
						ktree<k> largetree(large.data(), large.size());

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_hierarchical<k>(smalltree, largetree,
									out);
						}
						times["hierarchical"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_sequential<k>(smalltree, largetree, out);
						}
						times["sequential"] += timer.split();
//
//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							intersect_sorted<k>(small.data(), small.size(),
//									large.data(), large.size(), out);
//						}
//						times["sorted"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_hierarchical_bounded<k>(smalltree,
									largetree, out);
						}
						times["hierarchical_bounded"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_hierarchical_bounded_withoutLCA<k>(
									smalltree, largetree, out);
						}
						times["hierarchical_bounded_noLCA"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_sequential_bounded<k>(smalltree,
									largetree, out);
						}
						times["sequential_bounded"] += timer.split();

//						timer.reset();
//						for (uint32_t howmany = 0; howmany < REPETITION;
//								++howmany) {
//							intersect_sequential_bounded_withoutET(large.data(),
//									large.size(), small.data(), small.size(),
//									out);
//						}
//						times["sequential_bounded_noET"] += timer.split();

						// verification
//						it = multiset.begin();
//						vector<uint32_t> final_intersection = intersect(*it++,
//								*it++);
//						for (; it != multiset.end();)
//							final_intersection = intersect(final_intersection,
//									*it++);
//						if (out.size() != final_intersection.size()) {
//							std::cerr << "bad result!  " << std::endl;
//							return;
//						} else
//							printf("good!  ");
					} // for CASES
					for (auto time : times) {
						if (time.second != 0) {
							printf("%s: \e[31m%6.0f\e[0m  ", time.first.c_str(),
									(double) time.second / REPETITION / CASES);
						} // if !=0
					} // for times
					printf("\n");
					fflush(stdout);
				} // for num
			} // for size-ratio
		} // for intersection-ratio
	} // for minlength
}

template<karytree::indextype k = 3>
void run_hierarchical() {
	using namespace karytree;

	vector<uint32_t> out;

	uint32_t MaxBit = 31; // largest bit-length of element
	uint32_t minlength;
	size_t REPETITION = 1;
	size_t CASES = 20;

	WallClockTimer timer;
	vector<float> intersectionsratios = { 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
			0.70, 0.80, 0.90, 1.00 };
	vector<uint32_t> sizeratios = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
			2048, 4096, 8192 };

	std::map<std::string, size_t> times;
	for (uint32_t msb = 10; msb <= 14; ++msb) {
		minlength = 1U << msb;

		for (float ir : intersectionsratios) {
			printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
			for (uint32_t sr : sizeratios) {
				printf("  size ratio: \e[32m%4d\e[0m\n", sr);

				if (sr > 1000)
					REPETITION = 100;
				else
					REPETITION = 200;

				for (uint32_t num = 2; num < 3; num++) {
					times["hierarchical"] = 0;
					times["hierarchical_bounded"] = 0;
					times["hierarchical_bounded_new"] = 0;

					time_t t = time(nullptr);
					tm* _tm = localtime(&t);
					printf("%02d:%02d:%02d> num: \e[32m%2d\e[0m  \n",
							_tm->tm_hour, _tm->tm_min, _tm->tm_sec, num);

					for (uint32_t instance = 0; instance < CASES; instance++) {
						mySet multiset;

						ClusteredDataGenerator cdg;
						multiset = genMultipleSets(cdg, minlength, num,
								1U << MaxBit, static_cast<float>(sr), ir);

						// build k-ary tree
						auto it = multiset.begin();
						std::vector<keytype> small(std::move(*it++));
						std::vector<keytype> large(std::move(*it));

						ktree<3> smalltree(small.data(), small.size());
						ktree<3> largetree(large.data(), large.size());

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_hierarchical<k>(smalltree, largetree,
									out);
						}
						times["hierarchical"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_hierarchical_bounded<k>(smalltree,
									largetree, out);
						}
						times["hierarchical_bounded"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_hierarchical_bounded_new<k>(smalltree,
									largetree, out);
						}
						times["hierarchical_bounded_new"] += timer.split();

						// verification
//						it = multiset.begin();
//						vector<uint32_t> final_intersection = intersect(*it++,
//								*it++);
//						for (; it != multiset.end();)
//							final_intersection = intersect(final_intersection,
//									*it++);
//						if (out.size() != final_intersection.size()) {
//							std::cerr << "bad result!  " << std::endl;
//							return;
//						} else
//							printf("good!  ");
					} // for CASES
					for (auto time : times) {
						if (time.second != 0) {
							printf("%s: \e[31m%6.0f\e[0m  ", time.first.c_str(),
									(double) time.second / REPETITION / CASES);
						} // if !=0
					} // for times
					printf("\n");
					fflush(stdout);
				} // for num
			} // for size-ratio
		} // for intersection-ratio
	} // for minlength
}

template<karytree::indextype k = 3>
void test() {
//	using namespace karytree_trial;
	using namespace karytree;
	vector<uint32_t> out, out1;

	ClusteredDataGenerator cdg;
	uint32_t intersize = 10;

	for (int s = 8; s < 2560; s += 2) {
		for (int i = 0; i < 10; i++) {
			vector<uint32_t> inter = cdg.generate(intersize, 26);

			vector<uint32_t> small = unite(
					cdg.generate(static_cast<uint32_t>(s * 8 - inter.size()),
							s * 8), inter);
			small.resize(s);

			vector<uint32_t> large = unite(
					cdg.generate(static_cast<uint32_t>(s * 8 - inter.size()),
							s * 8), inter);
			large.resize(s);
//			small= {0, 1, 2, 3, 4, 5, 6,7/*, 8, 9, 10, 11, 12, 13, 14, 15, 16,17, 18, 19, 20, 21*/};
//			large= {0, 1, 2, 3, 4, 5, 6,7/*, 8, 9, 10, 11, 12, 13, 14, 15, 16,17, 18, 19, 20, 21*/};
			ktree<3> smalltree(small.data(), small.size());
			ktree<3> largetree(large.data(), large.size());
			vector<uint32_t> final_intersection = intersect(small, large);

			/* reorder elements in vectors*/
//			construct_tree_fast(small.data(), small.size());
//			construct_tree_fast(large.data(), large.size());
			//	traverse_tree(small.data(), 0, 0, small.size());
			//	traverse_tree(large.data(), 0, 0, large.size());
//			karytree::intersect_hierarchical_bounded<k>(small.data(),
//					small.size(), large.data(), large.size(), out);
//			karytree::intersect_hierarchical_bounded_new<k>(smalltree,
//					largetree, out1);
			karytree::intersect_sequential_bounded<k>(smalltree, largetree,
					out);
//			karytree::intersect_sequential_bounded<k>(small.data(),
//					small.size(), large.data(), large.size(), out);

			if (/*out != out1*/
			/*&&*/out.size() != final_intersection.size()) {
				std::cout << "small: ";
				for (uint32_t j = 0; j < small.size(); j++) {
					std::cout << small[j] << ",";
				}
				std::cout << std::endl << "large:: ";
				for (uint32_t j = 0; j < large.size(); j++) {
					std::cout << large[j] << ",";
				}
				std::cout << std::endl;
			}
		}
	}
	std::cout << "finish!" << std::endl;
}

void tree_construct_test() {
	using namespace karytree;
	std::vector<keytype> large;
	for (int i = 0; i < ktree<3>::PERFECT_SIZES[9]; i++)
		large.push_back(i + 1);
	large.resize(20);
	karytree::search_sorted_array<3>(1, large.data(), 19);

//	ktree<3> tree(large.data(), large.size());

	// lca test
//	for (int ld = 0; ld < 9; ld++)
//		for (int lr = 0; lr < ktree<3>::K_s[ld]; lr++)
//			for (int rd = 0; rd < 9; rd++)
//				for (int rr = 0; rr < ktree<3>::K_s[rd]; rr++)
//					tree.get_LCA(ld, lr, rd, rr);
//	for (int i = 0; i < tree.size(); i++)
//		std::cout << tree[i] << std::endl;
}

void generateMySets() {
	uint32_t MaxBit = 31; // largest bit-length of element
	uint32_t minlength;
	size_t REPETITION = 1;
	size_t CASES = 20;

	vector<float> intersectionsratios = { 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
			0.70, 0.80, 0.90, 1.00 };
	vector<uint32_t> sizeratios = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
			2048, 4096, 8192 };

	for (uint32_t msb = 10; msb <= 14; ++msb) {
		minlength = 1U << msb;

		for (float ir : intersectionsratios) {
			printf("intersection ratio: \e[32m%3.0f%%\e[0m\n", ir * 100);
			for (uint32_t sr : sizeratios) {
				printf("  size ratio: \e[32m%4d\e[0m\n", sr);

				if (sr > 1000)
					REPETITION = 100;
				else
					REPETITION = 200;

				for (uint32_t num = 2; num < 3; num++) {
					for (uint32_t instance = 0; instance < CASES; instance++) {
						mySet multiset;

						ClusteredDataGenerator cdg;
						std::cout << "we are here!" << std::endl;
						return;
						cdg.generate(minlength, 1 << MaxBit);
						multiset = genMultipleSets(cdg, minlength, num,
								1U << MaxBit, static_cast<float>(sr), ir);

					} // for CASES
				} // for num
			} // for size-ratio
		} // for intersection-ratio
	} // for minlength
}

int main(void) {

//	test();
//	tree_construct_test();
	generateMySets();

//	run<3>();
//	run_hierarchical();
	return 0;
}
