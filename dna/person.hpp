#pragma once

#include "sequence_buffer.hpp"

namespace dna
{

template<typename T>
concept HelixStream = requires(T a) {
	{ a.seek(1000L) };
	{ a.read() }; // -> std::convertible_to<sequence_buffer<ByteBuffer>>;
	{ a.size() } -> std::convertible_to<std::size_t>;
};

template<typename T>
concept Person = requires(T a) {
	{ a.chromosome(1) };
	{ a.chromosomes() } -> std::convertible_to<std::size_t>;
};

}

