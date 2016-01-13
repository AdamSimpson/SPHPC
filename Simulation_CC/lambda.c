#define HOST_LAMBDA [=]
#define DEVICE_LAMBDA [=]__device__
#define LAMBDA constexpr ((1==1)?(HOST_LAMBDA):(DEVICE_LAMBDA))

int main(int argc, char **argv) {
  LAMBDA
}
