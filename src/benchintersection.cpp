/*
 * benchintersection.cpp
 *
 *  Created on: 2018Äê1ÔÂ15ÈÕ
 *      Author: John
 */

#include "k-ary_search_tree_trial.hpp"
#include "k-ary_search_tree.hpp"

template<karytree_trial::indextype k = 3>
void run() {
	using namespace karytree_trial;

	vector<uint32_t> out;

	uint32_t MaxBit = 31; // largest bit-length of element
	uint32_t minlength;
	size_t REPETITION;
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
					times["sequential"] = 0;
					times["sorted"] = 0;
					times["hierarchical_bounded"] = 0;
					times["sequential_bounded"] = 0;
					times["sequential_bounded_noET"] = 0;
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

						construct_tree_fast(small.data(), small.size());
						construct_tree_fast(large.data(), large.size());

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_sequential<k>(small.data(), small.size(),
									large.data(), large.size(), out);
						}
						times["sequential"] += timer.split();
						timer.reset();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_sorted<k>(small.data(), small.size(),
									large.data(), large.size(), out);
						}
						times["sorted"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_hierarchical_bounded<k>(small.data(),
									small.size(), large.data(), large.size(),
									out);
						}
						times["hierarchical_bounded"] += timer.split();

						timer.reset();
						for (uint32_t howmany = 0; howmany < REPETITION;
								++howmany) {
							intersect_sequential_bounded<3>(small.data(),
									small.size(), large.data(), large.size(),
									out);
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

//						// verification
//						it = multiset.begin();
//						vector<uint32_t> final_intersection = intersect(*it++,
//								*it++);
//						for (; it != multiset.end(); it++)
//							final_intersection = intersect(final_intersection,
//									*it);
//						if (out != final_intersection) {
//							for (auto it : multiset) {
//								for (auto iit : it) {
//									std::cout << iit << ",";
//								}
//								std::cout << std::endl;
//							}
//							for (std::vector<uint32_t>::iterator iout =
//									out.begin(), ifinal =
//									final_intersection.begin();
//									iout != out.end(); iout++, ifinal++) {
//								if (*iout != *ifinal) {
//									std::cout
//											<< (iout - out.begin())
//													/ sizeof(uint32_t) << ","
//											<< *iout << std::endl;
//									ifinal++;
//								}
//							}
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

template<karytree_trial::indextype k = 3>
void test() {
	using namespace karytree_trial;

	std::vector<keytype> out;

	std::vector<keytype> small;
	std::vector<keytype> large;

	for (int i = 0; i < 9; i++)
		small.push_back(i + 1);
	for (int i = 0; i < 26; i++)
		large.push_back(i + 1);

	karytree::ktree<3> tree1(large.data(), large.size());
	karytree::ktree<5> tree2(large.data(), large.size());
	karytree::ktree<9> tree3(large.data(), large.size());
	karytree::ktree<17> tree4(large.data(), large.size());

	/*********************************/
	/* verification */
	/*********************************/
	vector<uint32_t> final_intersection = intersect(small, large);

	/*********************************/
	/* intersect via k-ary search */
	/*********************************/

	construct_tree(small.data(), small.size());
	construct_tree_fast(large.data(), large.size());

	intersect_sequential_bounded_withoutET(small.data(), small.size(),
			large.data(), large.size(), out);

	intersect_hierarchical_bounded<k>(small.data(), small.size(), large.data(),
			large.size(), out);
}

int main(void) {
	using namespace karytree_trial;

//	test<3>();
	run<3>();
	return 0;
}
