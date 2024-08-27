int main(int argc, char* argv[]) {
    int count = 30;
    if(argc == 2) {
        count = std::stoi(argv[1]);
    }

    for(int i = 0; i <= count; ++i) {
        std::string in_file = std::string("init_examples/contour_") + std::to_string(i);
        std::string out_file = std::string("res_examples/contour_") + std::to_string(i); 
    // FILES TO READ AND WRITE
        std::ifstream in(in_file);
        std::ofstream out(out_file);

