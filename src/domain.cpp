/**
 * SDCSIM - strand displacement simulator on SIMD||DNA architecture
 * @file domain.cpp
 * @brief Contain implementation of domain class
 * @author Lukáš Plevač <xpleva07@vut.cz>
 */

#include "domain.h"

Domain::Domain(DOMAIN_DT value) {
    this->set(value);
}

Domain::Domain() {
    this->set(0);
}

void Domain::set(DOMAIN_DT value) {
    this->value = value;
} 

DOMAIN_DT Domain::get() {
    return this->value;
}