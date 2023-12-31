#include <iostream>
#include <gmp.h>

typedef struct {
    mpz_t* a;
    mpz_t* b;
    mpz_t* c;
} PythagoreanTriple;

void printTriple(PythagoreanTriple& trip) {
    std::cout << "(";
    mpz_out_str(stdout, 10, *trip.a);
    std::cout << ", ";
    mpz_out_str(stdout, 10, *trip.b);
    std::cout << ", ";
    mpz_out_str(stdout, 10, *trip.c);
    std::cout << ")";
}

void reduceTriple(PythagoreanTriple& trip) {
    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcd(gcd, *trip.a, *trip.b);
    mpz_gcd(gcd, gcd, *trip.c);

    // if gcd of triple is >1
    if (mpz_cmp_ui(gcd, 1) > 0) {
        mpz_divexact(*trip.a, *trip.a, gcd);
        mpz_divexact(*trip.b, *trip.b, gcd);
        mpz_divexact(*trip.c, *trip.c, gcd);
    }
}

bool isPrimitive(PythagoreanTriple& trip) {
    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcd(gcd, *trip.a, *trip.b);
    mpz_gcd(gcd, gcd, *trip.c);

    return mpz_cmp_ui(gcd, 1) > 0;
}

bool isValidParms(mpz_t& n, mpz_t& m) {
    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcd(gcd, n, m);
    return mpz_cmp(n, m) < 0 && mpz_cmp_ui(gcd, 1) == 0 && !mpz_congruent_2exp_p(n, m, 1);
}

void setTripleByParams(mpz_t& n, mpz_t& m, PythagoreanTriple& trip) {
    mpz_t temp;
    mpz_init(temp);

    mpz_mul(*trip.a, n, m); // a = nm
    mpz_mul_2exp(*trip.a, *trip.a, 1); // a *= 2

    mpz_pow_ui(*trip.b, m, 2); // b = m^2
    mpz_pow_ui(temp, n, 2); // temp = n^2
    mpz_sub(*trip.b, *trip.b, temp); // b -= temp

    mpz_pow_ui(*trip.c, m, 2); // b = m^2
    mpz_pow_ui(temp, n, 2); // temp = n^2
    mpz_add(*trip.c, *trip.c, temp); // b -= temp
}

int main (int argc, char **argv) {

    mpz_t n, m, a, b, c;
    mpz_inits(n, m, a, b, c, NULL);

    PythagoreanTriple temp = {&a, &b, &c};

    mpz_set_str(n, "1", 10);
    mpz_set_str(m, "1", 10);

    std::cout << "n, m: (a, b, c)" << std::endl;
    std::cout << "---------------" << std::endl;

    while (mpz_cmp_ui(n, 10) <= 0) {
        mpz_add_ui(m, n, 1);
        while (mpz_cmp_ui(m, 10) <= 0) {
            if (isValidParms(n, m)) {
                mpz_out_str(stdout, 10, n);
                std::cout << ", ";
                mpz_out_str(stdout, 10, m);
                std::cout << ": ";
                setTripleByParams(n, m, temp);
                reduceTriple(temp);
                printTriple(temp);
                std::cout << std::endl;
            }
            
            mpz_add_ui(m, m, 1);
        }
        mpz_add_ui(n, n, 1);
    }

    return 0;
}
