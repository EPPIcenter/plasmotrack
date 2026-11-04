//
// CborFile writes CBOR item by item that a CBOR reader reads back as the value meant, with numbers in their shortest
// head, doubles always as 64-bit floats, and a buffer that flushes in the middle of a large array.
//

#include "gtest/gtest.h"

#include "core/io/CborFile.h"

#include <nlohmann/json.hpp>

#include <unistd.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using transmission_nets::core::io::CborFile;
namespace fs = std::filesystem;

namespace {
    struct TempFile {
        fs::path path = fs::temp_directory_path() / ("plasmotrack-cbor-" + std::to_string(::getpid()) + "-" + ::testing::UnitTest::GetInstance()->current_test_info()->name() + ".cbor");
        ~TempFile() { fs::remove(path); }
    };

    nlohmann::json readBack(const fs::path& path) {
        std::ifstream in(path, std::ios::binary);
        return nlohmann::json::from_cbor(in);
    }
}// namespace

TEST(CborFileTest, AMapOfArraysTextsNumbersAndBytesReadsBack) {
    TempFile file;
    {
        CborFile out(file.path);
        out.map(6);
        out.text("texts");
        out.texts({"a", "bc", std::string(300, 'x')});
        out.text("unsigned");
        out.numbers(std::vector<std::uint64_t>{0, 23, 24, 255, 256, 65535, 65536, 4294967295u, 4294967296u});
        out.text("signed");
        out.numbers(std::vector<std::int64_t>{-1, -24, -25, -256, -257, std::numeric_limits<std::int64_t>::min(), 7});
        out.text("doubles");
        out.numbers(std::vector<double>{0.0, 1.5, -2.25, 1e300});
        out.text("bytes");
        out.bytes({0x00, 0xff, 0x10});
        out.text("small");
        out.value({{"k", {1, 2, 3}}, {"s", "v"}});
        out.close();
    }

    const auto back = readBack(file.path);
    EXPECT_EQ(back.at("texts"), (std::vector<std::string>{"a", "bc", std::string(300, 'x')}));
    EXPECT_EQ(back.at("unsigned"), (std::vector<std::uint64_t>{0, 23, 24, 255, 256, 65535, 65536, 4294967295u, 4294967296u}));
    EXPECT_EQ(back.at("signed"), (std::vector<std::int64_t>{-1, -24, -25, -256, -257, std::numeric_limits<std::int64_t>::min(), 7}));
    EXPECT_EQ(back.at("doubles"), (std::vector<double>{0.0, 1.5, -2.25, 1e300}));
    const auto& bytes = back.at("bytes").get_binary();
    EXPECT_EQ(std::vector<std::uint8_t>(bytes.begin(), bytes.end()), (std::vector<std::uint8_t>{0x00, 0xff, 0x10}));
    EXPECT_EQ(back.at("small"), (nlohmann::json{{"k", {1, 2, 3}}, {"s", "v"}}));
}

// Doubles are written as 64-bit floats whatever their value, so a reader gets them back exactly.
TEST(CborFileTest, DoublesAreWrittenAsSixtyFourBitFloats) {
    TempFile file;
    {
        CborFile out(file.path);
        out.array(2);
        out.number(0.1);
        out.number(1.0);
        out.close();
    }
    std::ifstream in(file.path, std::ios::binary);
    const std::vector<std::uint8_t> bytes{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    EXPECT_EQ(bytes.size(), 1u + 2 * 9);
    EXPECT_EQ(bytes[1], 0xfb);
    EXPECT_EQ(bytes[10], 0xfb);
    EXPECT_EQ(readBack(file.path).at(0).get<double>(), 0.1);
}

// An array larger than the buffer flushes as it is written and still reads back whole.
TEST(CborFileTest, AnArrayLargerThanTheBufferReadsBack) {
    TempFile file;
    const std::vector<double> many(300000, 2.5);// 2.7 MB of doubles, past the 1 MB buffer
    {
        CborFile out(file.path);
        out.numbers(many);
        out.close();
    }
    const auto back = readBack(file.path);
    ASSERT_EQ(back.size(), many.size());
    EXPECT_EQ(back.front(), 2.5);
    EXPECT_EQ(back.back(), 2.5);
}
