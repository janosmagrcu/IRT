#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

template<std::size_t dimension>
class Faraday
{
    public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension> const& B, VecField<dimension>& Bnew)
    {
        auto const dx = m_grid->cell_size(Direction::X);
        if constexpr (dimension == 1)
        {
            // Bx is primal in x
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                auto const& Bx = B.x;
                
                auto& Bnewx = Bnew.x;
                
                Bnewx(ix) = Bx(ix);
            }
            
            // By, Bz are dual in x
            for (auto ix = m_grid->dual_dom_start(Direction::X);
                 ix <= m_grid->dual_dom_end(Direction::X); ++ix)
            {
                auto const& By = B.y;
                auto const& Bz = B.z;
                
                auto const& Ey = E.y;
                auto const& Ez = E.z;
                
                auto& Bnewy = Bnew.y;
                auto& Bnewz = Bnew.z;
                
                Bnewy(ix) = By(ix) - (Ez(ix)-Ez(ix-1))/(dx);
                Bnewz(ix) = Bz(ix) + (Ey(ix)-Ey(ix-1))/(dx);             
            }
        }
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_FARADAY_HPP
