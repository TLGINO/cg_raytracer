# g++ -fsanitize=address main.cpp && ./a.out && open result.ppm
# g++ main.cpp && ./a.out && open result.ppm
clear
g++ main.cpp --std=c++11 -O3 && ./a.out && open result.ppm