#include "catch.hpp"
#include "helix_utilities.hpp"
#include <sequence_buffer.hpp>

TEST_CASE("Sequence buffer compare all equal", "[helix utils]")
{
	const std::array<std::byte, 2> data = {
			dna::pack(dna::G, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};
	dna::sequence_buffer buf1(data), buf2(data);
	
	const auto mismatched_intervals_0 = helix::compare(buf1, buf2);
	const auto mismatched_intervals_1 = helix::compare(buf1, buf2, 8);

	REQUIRE(mismatched_intervals_0.size() == 0);
	REQUIRE(mismatched_intervals_1.size() == 0);
}

TEST_CASE("Sequence buffer compare 1 mismatch at beginning", "[helix utils]")
{
	const std::array<std::byte, 2> data1 = {
			dna::pack(dna::G, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals_0 = helix::compare(buf1, buf2);
	const auto mismatched_intervals_1 = helix::compare(buf1, buf2, 8);

	REQUIRE(mismatched_intervals_0.size() == 1);
	REQUIRE(mismatched_intervals_0[0].first == 0);
	REQUIRE(mismatched_intervals_0[0].second == 0);
	REQUIRE(mismatched_intervals_1.size() == 1);
	REQUIRE(mismatched_intervals_1[0].first == 8);
	REQUIRE(mismatched_intervals_1[0].second == 8);
}

TEST_CASE("Sequence buffer compare 1 mismatch with length 2 in middle", "[helix utils]")
{
	const std::array<std::byte, 2> data1 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::G, dna::G),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals_0 = helix::compare(buf1, buf2);
	const auto mismatched_intervals_1 = helix::compare(buf1, buf2, 8);

	REQUIRE(mismatched_intervals_0.size() == 1);
	REQUIRE(mismatched_intervals_0[0].first == 2);
	REQUIRE(mismatched_intervals_0[0].second == 3);
	REQUIRE(mismatched_intervals_1.size() == 1);
	REQUIRE(mismatched_intervals_1[0].first == 10);
	REQUIRE(mismatched_intervals_1[0].second == 11);
}

TEST_CASE("Sequence buffer compare 1 mismatch at end", "[helix utils]")
{
	const std::array<std::byte, 2> data1 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::G),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals_0 = helix::compare(buf1, buf2);
	const auto mismatched_intervals_1 = helix::compare(buf1, buf2, 8);

	REQUIRE(mismatched_intervals_0.size() == 1);
	REQUIRE(mismatched_intervals_0[0].first == 7);
	REQUIRE(mismatched_intervals_0[0].second == 7);
	REQUIRE(mismatched_intervals_1.size() == 1);
	REQUIRE(mismatched_intervals_1[0].first == 15);
	REQUIRE(mismatched_intervals_1[0].second == 15);
}

TEST_CASE("Sequence buffer compare all equal but different lengths", "[helix utils]")
{
	const std::vector<std::byte> data1 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
			dna::pack(dna::A, dna::A, dna::A, dna::A),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals_0 = helix::compare(buf1, buf2);
	const auto mismatched_intervals_1 = helix::compare(buf1, buf2, 8);

	REQUIRE(mismatched_intervals_0.size() == 1);
	REQUIRE(mismatched_intervals_0[0].first == 8);
	REQUIRE(mismatched_intervals_0[0].second == 11);
	REQUIRE(mismatched_intervals_1.size() == 1);
	REQUIRE(mismatched_intervals_1[0].first == 16);
	REQUIRE(mismatched_intervals_1[0].second == 19);
}

TEST_CASE("Sequence buffer compare 1 mismatch at end and different lengths", "[helix utils]")
{
	const std::vector<std::byte> data1 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::G),
			dna::pack(dna::A, dna::A, dna::G, dna::G),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals_0 = helix::compare(buf1, buf2);
	const auto mismatched_intervals_1 = helix::compare(buf1, buf2, 8);

	REQUIRE(mismatched_intervals_0.size() == 1);
	REQUIRE(mismatched_intervals_0[0].first == 7);
	REQUIRE(mismatched_intervals_0[0].second == 11);
	REQUIRE(mismatched_intervals_1.size() == 1);
	REQUIRE(mismatched_intervals_1[0].first == 15);
	REQUIRE(mismatched_intervals_1[0].second == 19);
}

TEST_CASE("Sequence buffer combine non-overlapping intervals", "[helix utils]")
{
    const std::vector<helix::interval_list> intervals = {
        {{1,1},{3,3}},
        {{2,2},{4,4}},
    };

    const auto mismatched_intervals = helix::combine(intervals);

    REQUIRE(mismatched_intervals.size() == 4);
    REQUIRE(mismatched_intervals[0].first == 1);
    REQUIRE(mismatched_intervals[0].second == 1);
    REQUIRE(mismatched_intervals[1].first == 2);
    REQUIRE(mismatched_intervals[1].second == 2);
    REQUIRE(mismatched_intervals[2].first == 3);
    REQUIRE(mismatched_intervals[2].second == 3);
    REQUIRE(mismatched_intervals[3].first == 4);
    REQUIRE(mismatched_intervals[3].second == 4);
}

TEST_CASE("Sequence buffer combine overlapping intervals", "[helix utils]")
{
    const std::vector<helix::interval_list> intervals = {
        {{1,2},{3,3}},
        {{3,5},{7,9}},
        {{9,12},{14,14}},
    };

    const auto mismatched_intervals = helix::combine(intervals);

    REQUIRE(mismatched_intervals.size() == 4);
    REQUIRE(mismatched_intervals[0].first == 1);
    REQUIRE(mismatched_intervals[0].second == 2);
    REQUIRE(mismatched_intervals[1].first == 3);
    REQUIRE(mismatched_intervals[1].second == 5);
    REQUIRE(mismatched_intervals[2].first == 7);
    REQUIRE(mismatched_intervals[2].second == 12);
    REQUIRE(mismatched_intervals[3].first == 14);
    REQUIRE(mismatched_intervals[3].second == 14);
}
