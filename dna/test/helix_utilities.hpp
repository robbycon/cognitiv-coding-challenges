#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>
#include <person.hpp>
#include <sequence_buffer.hpp>

namespace helix
{

using interval = std::pair<std::size_t,std::size_t>;
using interval_list = std::vector<interval>;

// This function returns an ascending list of [start_idx, end_idx) intervals (inclusive start, exclusive end)
// where differences in data occur between parameters a and b. The optional 'offset' parameter indicates where
// this sequence starts in the larger dataset (if applicable).
// Time Complexity: O(min(m, n)) where m and n are the sizes of sequences a and b.
// Space Complexity: O(k) where k is the number of mismatched intervals.
template<typename T>
//template<std::ranges::random_access_range T>
//	requires std::ranges::sized_range<T>
interval_list compare(const T& a, const T& b, const std::size_t offset = 0) {
	if (offset < 0)
		throw std::invalid_argument("offset cannot be less than 0");

	const std::size_t m = a.size(), n = b.size(), sz = std::min(m, n);
	interval_list mismatched_intervals;

	for (std::size_t i = 0; i < sz; ++i) {
		if (a[i] != b[i]) {
			const std::size_t start_idx = i;
			while (++i < sz && a[i] != b[i]);
			mismatched_intervals.emplace_back(start_idx + offset, i + offset);
		}
	}

	if (const std::size_t extra = std::max(m,n); extra > sz) {
		// If the last interval went to the end of the sequence, extend it
		if (!mismatched_intervals.empty() && (mismatched_intervals.back().second - offset) == sz)
			mismatched_intervals.back().second = extra + offset;
		// Otherwise, add the extra buffer as a new interval
		else
			mismatched_intervals.emplace_back(sz + offset, extra + offset);
	}

	return mismatched_intervals;
}

// This function takes a group of sorted interval_list objects and combines them into a single interval_list.
// A typical use case would be to call the helix::compare function over different segments of a larger set of
// comparison data. These results could have a case where one segment's final interval was [x, y), and the
// next adjacent segment's first interval was [y, z). The end result of this interval should actually be [x,z)
// instead of {[x, y), [y, z)}. This function handles these cases, and performs the operation with O(nlogk) time
// and O(m + k) space.
// Time Complexity: O(nlogk) where n is the total number of intervals and k is the number of interval lists.
// Space Complexity: O(m + k) where m is the total number of intervals returned to the caller (1 <= m <= n) and
// k is the number of interval lists.
interval_list combine(const std::vector<interval_list>& mismatched_intervals) {
	// Step 1: initialize a min heap to help combine different intervals from the different lists
	using pq_item = std::tuple<interval,int,std::size_t>; // this contains [the interval, the index of its parent list, the index within that list]
	const int k = mismatched_intervals.size();
	std::vector<pq_item> init;
	for (int i = 0; i < k; ++i) {
		if (!mismatched_intervals[i].empty())
			init.emplace_back(mismatched_intervals[i][0], i, 0);
	}

	const auto cmp = [](const pq_item& a, const pq_item& b) {
		return std::get<0>(a).first > std::get<0>(b).first;
	};

	std::priority_queue<pq_item,std::vector<pq_item>,decltype(cmp)> pq(init.begin(), init.end(), cmp);

	// Step 2: extract next mismatched interval, and combine with the previously seen one if applicable
	interval_list result;
	while (!pq.empty()) {
		const auto [interval, parent_idx, list_idx] = pq.top(); pq.pop();
		if (!result.empty() && result.back().second >= interval.first)
			result.back().second = std::max(result.back().second, interval.second);
		else
			result.emplace_back(std::move(interval));
		
		// Add the next interval from the interval_list of the current popped value to the min heap
		if (list_idx + 1 < mismatched_intervals[parent_idx].size())
			pq.emplace(mismatched_intervals[parent_idx][list_idx + 1], parent_idx, list_idx + 1);
	}

	return result;
}

// This function reads the entire chosen chromosome stream from the person, and pipes it to the ostringstream
// parameter. Parameter 'chromosome_idx' is zero-indexed.
template<dna::Person P>
void read(P& person, const std::size_t chromosome_idx, std::ostringstream& writer) {
	if (chromosome_idx < 0)
        throw std::invalid_argument("chromosome index cannot be negative");
    if (chromosome_idx >= person.chromosomes())
        throw std::invalid_argument("chromosome index specified does not exist in person");

	auto chromosome = person.chromosome(chromosome_idx);
	while (true) {
		auto buffer = chromosome.read();
		if (buffer.size() == 0) break;
		writer << buffer;
	}
}

// This function splits the string_view parameter 'sv' into segments of the specified 'window_size'. If
// 'window_size' is <= 0, then the return will be a vector of size 1 containing the entire 'sv' range.
// Otherwise, the vector will have 'window_size'-length ranges for all the elements except potentially
// the final one, which could be smaller if 'sv' is not divisible by 'window_size'.
std::vector<std::string_view> split(const std::string_view& sv, int window_size) {
	const int n = sv.size();
	// If the requested window size is negative or 0, create a single segment.
	if (window_size <= 0) window_size = n;
	std::vector<std::string_view> segments;
	for (int i = 0; i < n; i += window_size) {
		segments.push_back(sv.substr(i, window_size));
	}
	return segments;
}

// This function compares a specified chromosome of two people and returns a combined interval_list of
// all the mismatches. It first reads the entire stream of data, strips the telomeres (TBD), creates and
// dispatches 'window_size' segments to be compared, and combines the results to return to the caller.
template<dna::Person P>
interval_list compare_chromosome(const P& a, const P& b, const std::size_t chromosome_idx, int window_size = -1) {
	if (chromosome_idx < 0)
        throw std::invalid_argument("chromosome index cannot be negative");
    if (chromosome_idx >= a.chromosomes())
        throw std::invalid_argument("chromosome index specified does not exist in Person a");
    if (chromosome_idx >= b.chromosomes())
        throw std::invalid_argument("chromosome index specified does not exist in Person b");

	// Step 1: Read the chromosome streams from Persons 'a' and 'b'.
	const auto chrom_data_a = std::make_unique<std::ostringstream>(),
		chrom_data_b = std::make_unique<std::ostringstream>();
	read(a, chromosome_idx, *chrom_data_a); read(b, chromosome_idx, *chrom_data_b);

	// Step 2: Strip the telomeres from the beginning and end of the chromosomes.
	// Implementation TBD.

	// Step 3: Split the valid chromosome data into 'window_size' chunks, which could be sent to
	// their own thread or separate server for independent processing in step 4.
	const auto segments_a = split(chrom_data_a->view(), window_size),
		segments_b = split(chrom_data_b->view(), window_size);

	// Step 4: This loop would be replaced by a dispatcher to send these segments to their own
	// thread or server, which would call a single 'compare(segment_a, segment_b, curr_offset)'
	// and return their results to be collected here in 'mismatched_intervals'.
	const auto m = segments_a.size(), n = segments_b.size(), sz = std::max(m, n);
	std::vector<interval_list> mismatched_intervals;
	for (int i = 0; i < sz; ++i) {
		mismatched_intervals.emplace_back(
			compare(
				i < m ? segments_a[i] : std::string_view(nullptr, 0),
				i < n ? segments_b[i] : std::string_view(nullptr, 0),
				i * window_size
				)
			);
	}

	// Step 5: This combines the mismatched chromosome ranges from the separate threads/servers
	// in step 4 to return a unified result to the caller.
	return combine(mismatched_intervals);
}

} // namespace helix
