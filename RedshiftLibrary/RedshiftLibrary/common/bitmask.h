#ifndef _REDSHIFT_COMMON_BITMASK_
#define _REDSHIFT_COMMON_BITMASK_

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace NSEpic
{

template <class Flag, typename = typename std::enable_if<std::is_enum<Flag>::value>::type>

class CBitMask {

    // The type of the bitmask is the same as the enum's underlying type.
    using underlying_type = typename std::underlying_type<Flag>::type;

    // This method sets enum values.
    static constexpr underlying_type mask_value(Flag f) {
        return static_cast<underlying_type>(f); 
    }

public:

    // Default constructors create a bitmask with no flag selected.
    constexpr CBitMask() noexcept = default;

    constexpr CBitMask(std::nullptr_t) noexcept : m_mask(0) {}

    // Constructor creates a bitmask with just one flag set.
    constexpr CBitMask(Flag f) noexcept : m_mask(mask_value(f)) {}

    // Copy constructor.
    constexpr CBitMask(const CBitMask& bm) noexcept : m_mask(bm.m_mask) {}

    // Get the value of the bitmask corresponding to the flag.
    constexpr underlying_type GetMaskValue(Flag f) const {
        return m_mask & mask_value(f);
    }

    // Get the value of the bit corresponding to the flag.
    constexpr underlying_type GetBitValue(Flag f) const {
        return mask_value(f);
    }

    // Get the global value of the bitmask.
    constexpr underlying_type GetBitMaskValue() const noexcept {
        return m_mask;
    }

    /* Flag accesses */

    // Test if any flag is set.
    bool TestAnyFlag() const noexcept {
        return m_mask != 0;
    }

    // Test if no flag is set.
    bool TestNoFlag() const noexcept {
        return m_mask == 0;
    }

    // Test if all flags are set.
    bool TestAllFlags() const noexcept {
        underlying_type mask = m_mask;
        if (mask == 0) 
            return false;
        while (mask > 0) { 
            if ((mask & 1) == 0) 
                return false;                
            mask >>= 1;
        } 
        return true;
    }

    // Test if the bit corresponding to the flag is set.
    bool TestFlag(Flag f) const {
        return (m_mask & mask_value(f)) == mask_value(f);
    }

    // Returns the number of flags in the bitmask that are set.
    std::size_t count() const noexcept {
        underlying_type mask = m_mask;
        std::size_t count = 0; 
        while (mask > 0) { 
            count += (mask & 1); 
            mask >>= 1;
        }
        return count; 
    }

    // Returns the number of flags in the bitmask.
    std::size_t size() const noexcept {
        underlying_type mask = -1;
        std::size_t size = 1;
        while (mask >>= 1) {
            size ++;
        }
        return size;
    }

    /* Flag operations */

    // Set the bit corresponding to the flag.
    void SetFlag(Flag f) {
        m_mask |= mask_value(f);
    }

    // Set all flags.
    void SetAllFlags() noexcept {
        m_mask = -1;
    }

    // Reset the bit corresponding to the flag.
    void ResetFlag(Flag f) {
        m_mask &= ~mask_value(f);
    }

    // Reset all flags.
    void ResetAllFlags() noexcept {
        m_mask = 0;
    }

    // Toggle the bit corresponding to the flag.
    void ToggleFlag(Flag f) {
        m_mask ^= mask_value(f);
    }

    // Toggle all flags.
    void ToggleAllFlags() noexcept {
        unsigned int length = log2(m_mask) + 1;
        m_mask ^= (1 << length) - 1;
    }

    // Returns a string that represents all flags raised in the bitmask.
    std::string GetFlagsRaised() const noexcept {
        std::string flag_raised;
        if (m_mask == 0) {
            flag_raised = "none";
        }
        else {
            unsigned int numf = 1;
            unsigned int flag = 1;
            while (m_mask >= flag) {
                if ((m_mask & flag) == flag)
                    flag_raised += std::to_string(numf); 
                numf ++;
                flag <<= 1;
            }
        }
        return flag_raised;
    }

private:

    underlying_type m_mask = 0;

};

// Creates a flag from multiple flags at the same time.
template <class Flag, typename = typename std::enable_if<std::is_enum<Flag>::value>::type>
inline constexpr Flag operator | (const Flag& l, const Flag& r) noexcept
{
    return static_cast<Flag>(static_cast<typename std::underlying_type<Flag>::type>(l) |
                             static_cast<typename std::underlying_type<Flag>::type>(r));
}


// Definition of five quality flags.
enum QFlags : std::uint32_t
{
    QFLAGS_NONE  = 0,
    QFLAGS_ONE   = (1 << 0),
    QFLAGS_TWO   = (1 << 1),
    QFLAGS_THREE = (1 << 2),
    QFLAGS_FOUR  = (1 << 3),
    QFLAGS_FIVE  = (1 << 4)
};
using QualityFlags = CBitMask<enum QFlags>;

}

#endif
