#include "vecfield.hpp"
#include "field.hpp"
#include "gridlayout.hpp"
#include "faraday.hpp"
#include "particle.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <vector>
#include <cmath>

void known_B_E()
{
    std::cout << "Running Faraday test...\n";
    
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};

    for (auto ix = layout->ghost_start(Quantity::By, Direction::X); ix <= layout->ghost_end(Quantity::By, Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, static_cast<int>(ix));
        B.y(ix) = 2.0; // known By
        B.z(ix) = -1.0; // known Bz
    }

    for (auto ix = layout->ghost_start(Quantity::Ey, Direction::X); ix <= layout->ghost_end(Quantity::Ey, Direction::X);
         ++ix)
    {
        auto x = layout->coordinate(Direction::X, Quantity::By, static_cast<int>(ix));
        E.y(ix) = std::cos(x); // known Ey
        E.z(ix) = std::sin(x); // known Ez
    }

    Faraday<dimension> farad_test{layout};
    farad_test(E, B, Bnew);  
    {
        std::string filename = "known_B_E.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/Ex", E.x.data());
        file.createDataSet("/Ey", E.y.data());
        file.createDataSet("/Ez", E.z.data());
        
        file.createDataSet("/Bx", B.x.data());
        file.createDataSet("/By", B.y.data());
        file.createDataSet("/Bz", B.z.data());
        
        file.createDataSet("/Bnewx", Bnew.x.data());
        file.createDataSet("/Bnewy", Bnew.y.data());
        file.createDataSet("/Bnewz", Bnew.z.data());
    }
}

int main()
{
    known_B_E();
}
