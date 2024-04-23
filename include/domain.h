/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file domain.h
 * @brief Contain headers for domain class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#pragma once

#include <cstdint>
#include <limits.h>

#define DOMAIN_DT uint16_t
#define IS_COMPLEMENTARY(X) (X > (USHRT_MAX >> 1))
#define NORMALIZE_DOMAIN(X) (IS_COMPLEMENTARY(X) ? ~X : X)

/**
 * @brief Class reprezenting Domin in strand displacement system 
 */
class Domain {
    public:

        /**
         * Class conscructor
         * @param value inicialization value of domain
         */
        Domain(DOMAIN_DT value);

        /**
         * Class conscructor without init value
         */
        Domain();

        /**
         * Set value of current domain
         * @param value value to set 
         */
        void set(DOMAIN_DT value);

        /**
         * Get lowlevel value of domain
         */
        DOMAIN_DT get();

        /**
         * Get complement of domain
         */
        Domain operator~() {
            return Domain(~this->value);
        }

        bool operator==(Domain d2) {
            return (d2.get() == this->value);
        }

        bool operator!=(Domain d2) {
            return (d2.get() != this->value);
        }

    private:
        DOMAIN_DT value;
};