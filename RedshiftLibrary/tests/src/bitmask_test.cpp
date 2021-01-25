#include <RedshiftLibrary/common/bitmask.h>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;


BOOST_AUTO_TEST_SUITE(test_bitmask)

enum class Colors : uint8_t {
                        //Decimal  //Binary
    Yellow = (1 << 0),  // 1       00000001
    Red    = (1 << 1),  // 2       00000010
    Green  = (1 << 2),  // 4       00000100
    Blue   = (1 << 3)   // 8       00001000
};

enum class RainbowColors : uint16_t {
    Red    = (1 << 0),
    Orange = (1 << 1),
    Yellow = (1 << 2),
    Green  = (1 << 3),
    Blue   = (1 << 4),
    Indigo = (1 << 5),
    Violet = (1 << 6)
};


BOOST_AUTO_TEST_CASE(BITMASK_COLORS){

CBitMask<Colors> flags;

BOOST_CHECK_EQUAL(sizeof(flags), sizeof(uint8_t));
BOOST_CHECK_EQUAL(alignof(flags), alignof(uint8_t));
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow), false);
BOOST_CHECK_EQUAL(flags.GetBitValue(Colors::Yellow), 1);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Red), false);
BOOST_CHECK_EQUAL(flags.GetBitValue(Colors::Red), 2);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Green), false);
BOOST_CHECK_EQUAL(flags.GetBitValue(Colors::Green), 4);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Blue), false);
BOOST_CHECK_EQUAL(flags.GetBitValue(Colors::Blue), 8);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), 0);

// Set blue and yellow flags
flags.SetFlag(Colors::Blue);
flags.SetFlag(Colors::Yellow);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Blue), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Green), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Red), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow), true);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Blue|Colors::Yellow));

// Toggle blue flag
flags.ToggleFlag(Colors::Blue);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Blue), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Green), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Red), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow), true);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Yellow));

// Toggle red flag
flags.ToggleFlag(Colors::Red);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Blue), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Green), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Red), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow), true);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Red|Colors::Yellow));

// Reset yellow flag
flags.ResetFlag(Colors::Yellow);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Blue), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Green), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Red), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow), false);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Red));

// Toggle green flag
flags.ToggleFlag(Colors::Green);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Blue), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Green), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Red), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow), false);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Green|Colors::Red));

// Reset/Toggle green and red flags
flags.ResetFlag(Colors::Green);
flags.ToggleFlag(Colors::Red);

BOOST_CHECK_EQUAL(flags.TestAnyFlag(), false);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow|Colors::Red|Colors::Green|Colors::Blue), false);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), 0);

// Set blue flag
flags.SetFlag(Colors::Blue);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Blue), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow|Colors::Red|Colors::Green), false);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Yellow), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Red), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Green), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Blue), flags.GetBitValue(Colors::Blue));
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Blue));

// Set green flag
flags.SetFlag(Colors::Green);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Green|Colors::Blue), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow|Colors::Red), false);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Yellow), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Red), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Green), flags.GetBitValue(Colors::Green));
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Blue), flags.GetBitValue(Colors::Blue));
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Green|Colors::Blue));

// Set red flag
flags.SetFlag(Colors::Red);

BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Red|Colors::Green|Colors::Blue), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow), false);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Yellow), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Red), flags.GetBitValue(Colors::Red));
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Green), flags.GetBitValue(Colors::Green));
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Blue), flags.GetBitValue(Colors::Blue));
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Red|Colors::Green|Colors::Blue));

// Set yellow flag
flags.SetFlag(Colors::Yellow);

BOOST_CHECK_EQUAL(flags.TestAllFlags(), true);
BOOST_CHECK_EQUAL(flags.TestFlag(Colors::Yellow|Colors::Red|Colors::Green|Colors::Blue), true);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Yellow), flags.GetBitValue(Colors::Yellow));
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Red), flags.GetBitValue(Colors::Red));
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Green), flags.GetBitValue(Colors::Green));
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Blue), flags.GetBitValue(Colors::Blue));
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(Colors::Yellow|Colors::Red|Colors::Green|Colors::Blue));
BOOST_CHECK_EQUAL(flags.GetBitValue(Colors::Yellow|Colors::Red|Colors::Green|Colors::Blue), flags.GetBitValue(Colors::Yellow)+flags.GetBitValue(Colors::Red)+flags.GetBitValue(Colors::Green)+flags.GetBitValue(Colors::Blue));

// Reset all flags
flags.ResetAllFlags();

BOOST_CHECK_EQUAL(flags.TestNoFlag(), true);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetMaskValue(Colors::Yellow|Colors::Red|Colors::Green|Colors::Blue));
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Yellow|Colors::Red|Colors::Green|Colors::Blue), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Yellow), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Red), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Green), 0);
BOOST_CHECK_EQUAL(flags.GetMaskValue(Colors::Blue), 0);
BOOST_CHECK_EQUAL(flags.count(), 0);
BOOST_CHECK_EQUAL(flags.size(), 8);

}

BOOST_AUTO_TEST_CASE(BITMASK_RAINBOWCOLORS){

CBitMask<RainbowColors> rflags{nullptr};

BOOST_CHECK_EQUAL(sizeof(rflags), sizeof(uint16_t));
BOOST_CHECK_EQUAL(alignof(rflags), alignof(uint16_t));
BOOST_CHECK_EQUAL(rflags.TestNoFlag(), true);
BOOST_CHECK_EQUAL(rflags.GetFlagsRaised(), std::string("none"));
BOOST_CHECK_EQUAL(rflags.GetBitMaskValue(), 0);
BOOST_CHECK_EQUAL(rflags.count(), 0);
BOOST_CHECK_EQUAL(rflags.size(), 16);

// Set violet, red, blue and indigo flags at the same time
rflags.SetFlag(RainbowColors::Violet | RainbowColors::Red | RainbowColors::Blue | RainbowColors::Indigo);

BOOST_CHECK_EQUAL(rflags.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Red), true); 
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Red|RainbowColors::Orange), false);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Orange), false);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Yellow), false);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Green), false);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Green|RainbowColors::Blue), false);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Blue), true);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Indigo), true);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Violet), true);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Red), 1 << 0);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Red|RainbowColors::Red), 1 << 0);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Red|RainbowColors::Orange), 1 << 0); 
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Orange), 0);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Yellow), 0);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Green), 0);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Green|RainbowColors::Green), 0);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Green|RainbowColors::Blue), 1 << 4);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Blue), 1 << 4);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Indigo), 1 << 5);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Violet), 1 << 6);
BOOST_CHECK_EQUAL(rflags.GetFlagsRaised(), std::string("1567"));
BOOST_CHECK_EQUAL(rflags.GetBitMaskValue(), 113);
BOOST_CHECK_EQUAL(rflags.count(), 4);
BOOST_CHECK_EQUAL(rflags.size(), 16);

// Toggle orange, yellow and green flags at the same time
rflags.ToggleFlag(RainbowColors::Orange | RainbowColors::Yellow | RainbowColors::Green);

BOOST_CHECK_EQUAL(rflags.TestAllFlags(), true);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Orange), true);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Yellow), true);
BOOST_CHECK_EQUAL(rflags.TestFlag(RainbowColors::Green), true);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Orange), 1 << 1);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Yellow), 1 << 2);
BOOST_CHECK_EQUAL(rflags.GetMaskValue(RainbowColors::Green), 1 << 3);
BOOST_CHECK_EQUAL(rflags.GetFlagsRaised(), std::string("1234567"));
BOOST_CHECK_EQUAL(rflags.GetBitMaskValue(), 127);
BOOST_CHECK_EQUAL(rflags.count(), 7);
BOOST_CHECK_EQUAL(rflags.size(), 16);


CBitMask<RainbowColors> rflags2(RainbowColors::Violet | RainbowColors::Green | RainbowColors::Yellow);

// Toggle red flag
rflags2.ToggleFlag(RainbowColors::Red);

BOOST_CHECK_EQUAL(sizeof(rflags2), sizeof(uint16_t));
BOOST_CHECK_EQUAL(alignof(rflags2), alignof(uint16_t));
BOOST_CHECK_EQUAL(rflags2.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(rflags2.TestAllFlags(), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Violet), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Green), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Yellow), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Red), true);
BOOST_CHECK_EQUAL(rflags2.GetFlagsRaised(), std::string("1347"));
BOOST_CHECK_EQUAL(rflags2.GetBitMaskValue(), 77);
BOOST_CHECK_EQUAL(rflags2.count(), 4);
BOOST_CHECK_EQUAL(rflags2.size(), 16);

// Set indigo and blue flags at the same time
rflags2.SetFlag(RainbowColors::Indigo | RainbowColors::Blue);

BOOST_CHECK_EQUAL(rflags2.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(rflags2.TestAllFlags(), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Indigo), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Blue), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Red), true);
BOOST_CHECK_EQUAL(rflags2.GetFlagsRaised(), std::string("134567"));
BOOST_CHECK_EQUAL(rflags2.GetBitMaskValue(), 125);
BOOST_CHECK_EQUAL(rflags2.count(), 6);
BOOST_CHECK_EQUAL(rflags2.size(), 16);

// Reset violet, indigo, blue and green flags at the same time
rflags2.ResetFlag(RainbowColors::Violet | RainbowColors::Indigo | RainbowColors::Blue | RainbowColors::Green);

BOOST_CHECK_EQUAL(rflags2.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(rflags2.TestAllFlags(), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Violet), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Indigo), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Blue), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Green), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Yellow), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Orange), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Red), true);
BOOST_CHECK_EQUAL(rflags2.GetFlagsRaised(), std::string("13"));
BOOST_CHECK_EQUAL(rflags2.GetBitMaskValue(), 5);
BOOST_CHECK_EQUAL(rflags2.count(), 2);
BOOST_CHECK_EQUAL(rflags2.size(), 16);

// Toggle violet, blue and yellow flags at the same time
rflags2.ToggleFlag(RainbowColors::Violet | RainbowColors::Blue | RainbowColors::Yellow);

BOOST_CHECK_EQUAL(rflags2.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(rflags2.TestAllFlags(), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Violet), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Indigo), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Blue), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Green), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Yellow), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Orange), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Red), true);
BOOST_CHECK_EQUAL(rflags2.GetFlagsRaised(), std::string("157"));
BOOST_CHECK_EQUAL(rflags2.GetBitMaskValue(), 81);
BOOST_CHECK_EQUAL(rflags2.count(), 3);
BOOST_CHECK_EQUAL(rflags2.size(), 16);

// Toggle all flags
rflags2.ToggleAllFlags();

BOOST_CHECK_EQUAL(rflags2.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(rflags2.TestAllFlags(), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Violet), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Indigo), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Blue), false);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Green), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Yellow), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Orange), true);
BOOST_CHECK_EQUAL(rflags2.TestFlag(RainbowColors::Red), false);
BOOST_CHECK_EQUAL(rflags2.GetFlagsRaised(), std::string("2346"));
BOOST_CHECK_EQUAL(rflags2.GetBitMaskValue(), 46);
BOOST_CHECK_EQUAL(rflags2.count(), 4);
BOOST_CHECK_EQUAL(rflags2.size(), 16);

// Reset all flags
rflags2.ResetAllFlags();

BOOST_CHECK_EQUAL(rflags2.TestNoFlag(), true);
BOOST_CHECK_EQUAL(rflags2.GetFlagsRaised(), std::string("none"));
BOOST_CHECK_EQUAL(rflags2.GetBitMaskValue(), 0);
BOOST_CHECK_EQUAL(rflags2.count(), 0);
BOOST_CHECK_EQUAL(rflags2.size(), 16);

// Set all flags
rflags2.SetAllFlags();

BOOST_CHECK_EQUAL(rflags2.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(rflags2.TestAllFlags(), true);
BOOST_CHECK_EQUAL(rflags2.GetFlagsRaised(), std::string("12345678910111213141516"));
BOOST_CHECK_EQUAL(rflags2.count(), 16);
BOOST_CHECK_EQUAL(rflags2.count(), rflags2.size());

}

BOOST_AUTO_TEST_CASE(QUALITY_FLAGS){ 

QualityFlags flags;

BOOST_CHECK_EQUAL(sizeof(flags), sizeof(uint32_t));
BOOST_CHECK_EQUAL(alignof(flags), alignof(uint32_t));
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_NONE), true);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_TWO), false);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_FIVE), false);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_ONE|QFLAGS_TWO|QFLAGS_THREE), false);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_THREE|QFLAGS_FOUR|QFLAGS_FIVE), false);
BOOST_CHECK_EQUAL(flags.GetFlagsRaised(), std::string("none"));
BOOST_CHECK_EQUAL(flags.TestAllFlags(), false);
BOOST_CHECK_EQUAL(flags.TestAnyFlag(), false);
BOOST_CHECK_EQUAL(flags.TestNoFlag(), true);
BOOST_CHECK_EQUAL(flags.count(), 0);
BOOST_CHECK_EQUAL(flags.size()-flags.count(), 32); 
 
// Set all quality flags
flags.SetAllFlags();

BOOST_CHECK_EQUAL(flags.TestNoFlag(), false);
BOOST_CHECK_EQUAL(flags.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(flags.TestAllFlags(), true);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), 4294967295);
BOOST_CHECK_EQUAL(flags.size(), 32);
BOOST_CHECK_EQUAL(flags.size()-flags.count(), 0);

// Toggle first quality flag
flags.ToggleFlag(QFLAGS_ONE);

BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_ONE), false);
BOOST_CHECK_EQUAL(flags.TestAllFlags(), false);
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), 4294967294);
BOOST_CHECK_EQUAL(flags.size(), 32);
BOOST_CHECK_EQUAL(flags.size()-flags.count(), 1);

// Reset all quality flags
flags.ResetAllFlags();

BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_NONE), true);
BOOST_CHECK_EQUAL(flags.GetFlagsRaised(), std::string("none"));
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), flags.GetBitValue(QFLAGS_NONE));
BOOST_CHECK_EQUAL(flags.size()-flags.count(), 32);

// Set multiple quality flags
flags.SetFlag(QFLAGS_TWO | QFLAGS_FOUR | QFLAGS_FIVE);

BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_ONE|QFLAGS_THREE), false);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_TWO|QFLAGS_FOUR|QFLAGS_FIVE), true);
BOOST_CHECK_EQUAL(flags.GetFlagsRaised(), std::string("245"));
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), 26);
BOOST_CHECK_EQUAL(flags.TestNoFlag(), false);
BOOST_CHECK_EQUAL(flags.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(flags.count(), 3);

// Set multiple quality flags
flags.SetFlag(QFLAGS_ONE | QFLAGS_THREE);

BOOST_CHECK_EQUAL(flags.GetFlagsRaised(), std::string("12345")); 
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), 31);
BOOST_CHECK_EQUAL(flags.TestAllFlags(), true);
BOOST_CHECK_EQUAL(flags.count(), 5);
BOOST_CHECK_EQUAL(flags.size(), 32);

// Copy constructor
QualityFlags qf(flags);

BOOST_CHECK_EQUAL(sizeof(qf), sizeof(uint32_t));
BOOST_CHECK_EQUAL(alignof(qf), alignof(uint32_t));
BOOST_CHECK_EQUAL(qf.GetFlagsRaised(), std::string("12345")); 
BOOST_CHECK_EQUAL(qf.GetBitMaskValue(), flags.GetBitMaskValue());
BOOST_CHECK_EQUAL(qf.count(), flags.count());
BOOST_CHECK_EQUAL(qf.size(), flags.size());

// Toggle multiple quality flags
qf.ToggleFlag(QFLAGS_TWO | QFLAGS_FOUR);

BOOST_CHECK_EQUAL(qf.TestAnyFlag(), true);
BOOST_CHECK_EQUAL(qf.TestFlag(QFLAGS_ONE), true);
BOOST_CHECK_EQUAL(qf.TestFlag(QFLAGS_TWO), false);
BOOST_CHECK_EQUAL(qf.TestFlag(QFLAGS_THREE), true);
BOOST_CHECK_EQUAL(qf.TestFlag(QFLAGS_FOUR), false);
BOOST_CHECK_EQUAL(qf.TestFlag(QFLAGS_FIVE), true);
BOOST_CHECK_EQUAL(qf.GetFlagsRaised(), std::string("135"));
BOOST_CHECK_EQUAL(qf.GetBitMaskValue(), flags.GetBitMaskValue()-qf.GetBitValue(QFLAGS_TWO|QFLAGS_FOUR));
BOOST_CHECK_EQUAL(qf.count()-flags.count(), -2);
BOOST_CHECK_EQUAL(qf.size()-flags.size(), 0);

// Assignment operation with multiple quality flags 
flags = static_cast<QualityFlags>(QFLAGS_THREE | QFLAGS_FIVE);

BOOST_CHECK_EQUAL(sizeof(flags), sizeof(uint32_t));
BOOST_CHECK_EQUAL(alignof(flags), alignof(uint32_t));
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_ONE), false);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_TWO), false);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_THREE), true);
BOOST_CHECK_EQUAL(flags.GetMaskValue(QFLAGS_THREE), flags.GetBitValue(QFLAGS_THREE));
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_FOUR), false);
BOOST_CHECK_EQUAL(flags.TestFlag(QFLAGS_FIVE), true);
BOOST_CHECK_EQUAL(flags.GetMaskValue(QFLAGS_FIVE), flags.GetBitValue(QFLAGS_FIVE));
BOOST_CHECK_EQUAL(flags.GetFlagsRaised(), std::string("35"));
BOOST_CHECK_EQUAL(flags.GetBitMaskValue(), 20);
BOOST_CHECK_EQUAL(flags.count(), 2);
BOOST_CHECK_EQUAL(flags.size(), 32);

}

BOOST_AUTO_TEST_SUITE_END() 
