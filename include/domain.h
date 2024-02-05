#pragma once

#include <cstdint>

#define DOMAIN_DT uint16_t

class Domain {
    public:
        Domain(DOMAIN_DT value);

        Domain();

        void set(DOMAIN_DT value);

        DOMAIN_DT get();

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