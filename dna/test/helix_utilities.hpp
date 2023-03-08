#pragma once

#include <algorithm>
#include <queue>
#include <utility>
#include <vector>
#include <sequence_buffer.hpp>

namespace helix
{

using interval = std::pair<uint64_t,uint64_t>;
using interval_list = std::vector<interval>;

// This function returns an ascending list of [start_idx, end_idx] intervals (inclusive) where differences
// in data occur between parameters a and b.
template<dna::ByteBuffer T>
interval_list compare(const dna::sequence_buffer<T>& a, const dna::sequence_buffer<T>& b) {
	const uint64_t m = a.size(), n = b.size(), sz = std::min(m, n);
	std::vector<std::pair<uint64_t,uint64_t>> mismatched_intervals;

	for (uint64_t i = 0; i < sz; ++i) {
		if (a[i] != b[i]) {
			const uint64_t start_idx = i;
			while (i < sz - 1 && a[i + 1] != b[i + 1]) ++i;
			mismatched_intervals.emplace_back(start_idx, i);
		}
	}

	if (const uint64_t extra = std::max(m,n); extra > sz) {
		// If the last interval went to the end of the sequence, extend it
		if (!mismatched_intervals.empty() && mismatched_intervals.back().second == sz - 1)
			mismatched_intervals.back().second = extra - 1;
		// Otherwise, add the extra buffer as a new interval
		else
			mismatched_intervals.emplace_back(sz, extra - 1);
	}

	return mismatched_intervals;
}

interval_list combine(const std::vector<interval_list>& intervals) {
	// Step 1: initialize a min heap to help combine different intervals from the different lists
	using pq_item = std::tuple<interval,int,int>; // this contains [the interval, the index of its parent list, the index within that list]
	const int n = intervals.size();
	std::vector<pq_item> init;
	for (int i = 0; i < n; ++i) {
		if (!intervals[i].empty())
			init.emplace_back(intervals[i][0], i, 0);
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
		
		if (list_idx + 1 < intervals[parent_idx].size())
			pq.emplace(intervals[parent_idx][list_idx + 1], parent_idx, list_idx + 1);
	}

	return result;
}

} // namespace helix
