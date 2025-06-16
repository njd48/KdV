

#include "KdV.cpp"

int main () {

    int N = 16;

    KdV_env sys(N);

    sys.dt      = 0.000001;
    sys.t_final = 0.06;
    sys.t_write = 0.03;

    std::cout << "sim begin.\n";

    sys.set_initial();
    sys.run();

    std::cout << "sim finished.\n";

    return 0;
}