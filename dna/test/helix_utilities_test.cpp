#include "catch.hpp"
#include "fake_person.hpp"
#include "helix_utilities.hpp"
#include <sequence_buffer.hpp>
#include <iostream>
#include <string_view>

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

std::vector<std::byte> to_byte(std::initializer_list<int> il) {
    const int n = il.size();
    std::vector<std::byte> data(n);
    int i = 0;
    for (auto it = il.begin(); it != il.end(); ++it, ++i)
        data[i] = static_cast<std::byte>(*it);
    return data;
}

TEST_CASE("Read Chromosome from Person", "[helix utils]")
{
    const auto data = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x64});
    const std::size_t chunk_size = 4;
    fake_person person(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size);

    std::ostringstream ss;
    helix::read(person, 0, ss);
    
    REQUIRE(ss.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");
}

TEST_CASE("Split Chromosome from Person negative window size", "[helix utils]")
{
    const auto data = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x64});
    const std::size_t chunk_size = 4;
    fake_person person(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size);

    std::ostringstream ss;
    helix::read(person, 0, ss);
    
    REQUIRE(ss.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");

    const auto segments = helix::split(ss.view(), -1);

    REQUIRE(segments.size() == 1);
    REQUIRE(segments[0] == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");
}

TEST_CASE("Split Chromosome from Person 0 window size", "[helix utils]")
{
    const auto data = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x64});
    const std::size_t chunk_size = 4;
    fake_person person(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size);

    std::ostringstream ss;
    helix::read(person, 0, ss);
    
    REQUIRE(ss.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");

    const auto segments = helix::split(ss.view(), 0);

    REQUIRE(segments.size() == 1);
    REQUIRE(segments[0] == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");
}

TEST_CASE("Split Chromosome from Person 4 segments", "[helix utils]")
{
    const auto data = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x64});
    const std::size_t chunk_size = 4;
    fake_person person(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size);

    std::ostringstream ss;
    helix::read(person, 0, ss);
    
    REQUIRE(ss.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");

    const auto segments = helix::split(ss.view(), 8);

    REQUIRE(segments.size() == 4);
    REQUIRE(segments[0] == "CCGGTGAT");
    REQUIRE(segments[1] == "ATTGATTT");
    REQUIRE(segments[2] == "GATCTGTC");
    REQUIRE(segments[3] == "CATCCGCA");
}

TEST_CASE("Read Chromosome from 2 Persons and compare, all equal", "[helix utils]")
{
    const auto data = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x64});
    const std::size_t chunk_size = 4;
    fake_person person1(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size),
    person2(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size);
    
    std::ostringstream ss1, ss2;
    helix::read(person1, 0, ss1);
    helix::read(person2, 0, ss2);
    
    REQUIRE(ss1.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");
    REQUIRE(ss2.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");

    const auto mismatched_intervals = helix::compare_chromosome(person1, person2, 0);

    REQUIRE(mismatched_intervals.size() == 0);
}

TEST_CASE("Read Chromosome from 2 Persons and compare, 1 mismatch at end", "[helix utils]")
{
    const auto data1 = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x64}),
    data2 = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x65});
    const std::size_t chunk_size = 4;
    fake_person person1(std::array<std::vector<std::byte>, 23> {
        data1, data1, data1, data1, data1, data1, data1, data1,
        data1, data1, data1, data1, data1, data1, data1, data1,
        data1, data1, data1, data1, data1, data1, data1
    }, chunk_size),
    person2(std::array<std::vector<std::byte>, 23> {
        data2, data2, data2, data2, data2, data2, data2, data2,
        data2, data2, data2, data2, data2, data2, data2, data2,
        data2, data2, data2, data2, data2, data2, data2
    }, chunk_size);
    
    std::ostringstream ss1, ss2;
    helix::read(person1, 0, ss1);
    helix::read(person2, 0, ss2);
    
    REQUIRE(ss1.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCA");
    REQUIRE(ss2.view() == "CCGGTGATATTGATTTGATCTGTCCATCCGCC");

    const auto mismatched_intervals = helix::compare_chromosome(person1, person2, 0);

    REQUIRE(mismatched_intervals.size() == 1);
    REQUIRE(mismatched_intervals[0].first == 31);
    REQUIRE(mismatched_intervals[0].second == 31);
}

TEST_CASE("Person compare all equal", "[helix utils]")
{
    const auto data = to_byte({0x5a, 0xe3, 0x3e, 0x3f, 0x8d, 0xed, 0x4d, 0x64});
    const std::size_t chunk_size = 4;
    const fake_person person1(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size),
    person2(std::array<std::vector<std::byte>, 23> {
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data, data,
        data, data, data, data, data, data, data
    }, chunk_size);

    const auto mismatched_intervals = helix::compare_chromosome(person1, person2, 0);

    REQUIRE(mismatched_intervals.size() == 0);
}
