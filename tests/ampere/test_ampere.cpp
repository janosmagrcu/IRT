#include "vecfield.hpp"
#include "field.hpp"
#include "gridlayout.hpp"
#include "ampere.hpp"
#include "particle.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <vector>
#include <cmath>

void sine_B()
{
    std::cout << "Running Ampere test...\n";
    
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};


    for (auto ix = layout->ghost_start(Quantity::By, Direction::X); ix <= layout->ghost_end(Quantity::By, Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, static_cast<int>(ix));
        B.y(ix) = std::cos(x); // known By
        B.z(ix) = std::sin(x); // known Bz
    }

    Ampere<dimension> amp_test{layout};
    amp_test(B, J);  
    {
        std::string filename = "sine_B.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/Bx", B.x.data());
        file.createDataSet("/By", B.y.data());
        file.createDataSet("/Bz", B.z.data());
        
        file.createDataSet("/Jx", J.x.data());
        file.createDataSet("/Jy", J.y.data());
        file.createDataSet("/Jz", J.z.data());
    }
}

int main()
{
    sine_B();
}
