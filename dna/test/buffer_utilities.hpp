#pragma once

#include <algorithm>
#include <utility>
#include <vector>
#include <sequence_buffer.hpp>

// This function returns an ascending list of [start_idx, end_idx] intervals (inclusive) where differences
// in data occur between parameters a and b.
template<dna::ByteBuffer T>
std::vector<std::pair<uint64_t,uint64_t>> compare(const dna::sequence_buffer<T>& a, const dna::sequence_buffer<T>& b) {
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
