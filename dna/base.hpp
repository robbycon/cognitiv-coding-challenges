#pragma once

#include <cstddef>

#include <ostream>

namespace dna
{

enum class base
{
	adenine,
	cytosine,
	guanine,
	thymine
};

static constexpr base A = base::adenine;
static constexpr base C = base::cytosine;
static constexpr base G = base::guanine;
static constexpr base T = base::thymine;

constexpr std::byte pack(base first, base second, base third, base fourth)
{
	return (static_cast<std::byte>(first) << 6) |
					(static_cast<std::byte>(second) << 4) |
					(static_cast<std::byte>(third) << 2) |
					(static_cast<std::byte>(fourth));
}

using packed_bases = std::array<base, 4>;
using packed_size = std::tuple_size<packed_bases>;

constexpr packed_bases unpack(std::byte b)
{
	return packed_bases({
			static_cast<dna::base>((b >> 6) & std::byte{0x3}),
			static_cast<dna::base>((b >> 4) & std::byte{0x3}),
			static_cast<dna::base>((b >> 2) & std::byte{0x3}),
			static_cast<dna::base>(b & std::byte{0x3}) });
}

inline std::ostream& operator<<(std::ostream& os, base v)
{
	switch (v)
	{
		case base::adenine:
			return os << 'A';
		case base::cytosine:
			return os << 'C';
		case base::guanine:
			return os << 'G';
		default:
			return os << 'T';
	}
}

}
