#include "catch.hpp"
#include <array>
#include "helix_utilities.hpp"
#include "sequence_buffer.hpp"

TEST_CASE("Can use a Sequence Buffer", "[seqbuf]")
{
	std::array<std::byte, 2> data = {
			dna::pack(dna::G, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};

	dna::sequence_buffer buf(data);
	REQUIRE(buf[0] == dna::G);
	REQUIRE(buf[1] == dna::A);
	REQUIRE(buf[2] == dna::C);
	REQUIRE(buf[3] == dna::T);
	REQUIRE(buf[4] == dna::A);
	REQUIRE(buf[5] == dna::A);
	REQUIRE(buf[6] == dna::G);
	REQUIRE(buf[7] == dna::C);
}

TEST_CASE("Can use an iterator", "[seqbuf]")
{
	std::array<std::byte, 2> data = {
			dna::pack(dna::G, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};

	dna::sequence_buffer buf(data);
	INFO("SEQUENCE: " << buf);

	std::vector<dna::base> bases;
	for (auto i : buf)
		bases.emplace_back(i);

	REQUIRE(bases[0] == dna::G);
	REQUIRE(bases[1] == dna::A);
	REQUIRE(bases[2] == dna::C);
	REQUIRE(bases[3] == dna::T);
	REQUIRE(bases[4] == dna::A);
	REQUIRE(bases[5] == dna::A);
	REQUIRE(bases[6] == dna::G);
	REQUIRE(bases[7] == dna::C);

}

TEST_CASE("Sequence buffer compare all equal", "[stream]")
{
	const std::array<std::byte, 2> data = {
			dna::pack(dna::G, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};
	dna::sequence_buffer buf1(data), buf2(data);
	
	const auto mismatched_intervals = helix::compare(buf1, buf2);
	REQUIRE(mismatched_intervals.size() == 0);
}

TEST_CASE("Sequence buffer compare 1 mismatch at beginning", "[stream]")
{
	const std::array<std::byte, 2> data1 = {
			dna::pack(dna::G, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals = helix::compare(buf1, buf2);
	REQUIRE(mismatched_intervals.size() == 1);
	REQUIRE(mismatched_intervals[0].first == 0);
	REQUIRE(mismatched_intervals[0].second == 0);
}

TEST_CASE("Sequence buffer compare 1 mismatch with length 2 in middle", "[stream]")
{
	const std::array<std::byte, 2> data1 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::G, dna::G),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals = helix::compare(buf1, buf2);
	REQUIRE(mismatched_intervals.size() == 1);
	REQUIRE(mismatched_intervals[0].first == 2);
	REQUIRE(mismatched_intervals[0].second == 3);
}

TEST_CASE("Sequence buffer compare 1 mismatch at end", "[stream]")
{
	const std::array<std::byte, 2> data1 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::C),
	}, data2 = {
			dna::pack(dna::T, dna::A, dna::C, dna::T),
			dna::pack(dna::A, dna::A, dna::G, dna::G),
	};
	dna::sequence_buffer buf1(data1), buf2(data2);
	
	const auto mismatched_intervals = helix::compare(buf1, buf2);
	REQUIRE(mismatched_intervals.size() == 1);
	REQUIRE(mismatched_intervals[0].first == 7);
	REQUIRE(mismatched_intervals[0].second == 7);
}

TEST_CASE("Sequence buffer compare all equal but different lengths", "[stream]")
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
	
	const auto mismatched_intervals = helix::compare(buf1, buf2);
	REQUIRE(mismatched_intervals.size() == 1);
	REQUIRE(mismatched_intervals[0].first == 8);
	REQUIRE(mismatched_intervals[0].second == 11);
}

TEST_CASE("Sequence buffer compare 1 mismatch at end and different lengths", "[stream]")
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
	
	const auto mismatched_intervals = helix::compare(buf1, buf2);
	REQUIRE(mismatched_intervals.size() == 1);
	REQUIRE(mismatched_intervals[0].first == 7);
	REQUIRE(mismatched_intervals[0].second == 11);
}
