#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include "vecfield.hpp"
#include "particle.hpp"

#include <cstddef>
#include <vector>


template<std::size_t dimension>
class Pusher
{
protected:
    std::shared_ptr<GridLayout<dimension>> layout_;
    double dt_;

public:
    Pusher(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : layout_(layout)
        , dt_(dt)
    {
    }

    virtual void operator()(std::vector<Particle<dimension>>& particles,
                            VecField<dimension> const& E, VecField<dimension> const& B)
        = 0;

    virtual ~Pusher() {}
};



template<std::size_t dimension>
class Boris : public Pusher<dimension>
{
public:
    Boris(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : Pusher<dimension>{layout, dt}
    {
    }


    void operator()(std::vector<Particle<dimension>>& particles, VecField<dimension> const& E,
                VecField<dimension> const& B) override
    {
        for (auto& particle : particles)
        {
            double const dx = this->layout_->cell_size(Direction::X);
            double const dt = this->dt_;
            
            particle.position[0] += 0.5 * particle.v[0] * dt;
                
            int const iCell = static_cast<int>(particle.position[0]/dx + this->layout_->dual_dom_start(Direction::X)); 
            double reminder = (particle.position[0] / dx) - iCell;

            double Ex = this->interpolate(E.x, iCell, reminder);
            double Ey = this->interpolate(E.y, iCell, reminder);
            double Ez = this->interpolate(E.z, iCell, reminder);
    
            double Bx = this->interpolate(B.x, iCell, reminder);
            double By = this->interpolate(B.y, iCell, reminder);
            double Bz = this->interpolate(B.z, iCell, reminder);
    
            double qmdt2 = (particle.charge / particle.mass) * 0.5 * dt;

            double vminus_x = particle.v[0] + qmdt2 * Ex;
            double vminus_y = particle.v[1] + qmdt2 * Ey;
            double vminus_z = particle.v[2] + qmdt2 * Ez;
    
            double tx = qmdt2 * Bx;
            double ty = qmdt2 * By;
            double tz = qmdt2 * Bz;
    
            double t2 = tx*tx + ty*ty + tz*tz;

            // painful vector manipulations

            double vprime_x = vminus_x + (vminus_y*tz - vminus_z*ty);
            double vprime_y = vminus_y + (vminus_z*tx - vminus_x*tz);
            double vprime_z = vminus_z + (vminus_x*ty - vminus_y*tx);
 
            double sx = 2.0 * tx / (1.0 + t2);
            double sy = 2.0 * ty / (1.0 + t2);
            double sz = 2.0 * tz / (1.0 + t2);

            double vplus_x = vminus_x + (vprime_y*sz - vprime_z*sy);
            double vplus_y = vminus_y + (vprime_z*sx - vprime_x*sz);
            double vplus_z = vminus_z + (vprime_x*sy - vprime_y*sx);
    
            particle.v[0] = vplus_x + qmdt2 * Ex;
            particle.v[1] = vplus_y + qmdt2 * Ey;
            particle.v[2] = vplus_z + qmdt2 * Ez;

            particle.position[0] += 0.5 * particle.v[0] * dt;
        }
    }

private:
    double interpolate(Field<dimension> const& field, int iCell, double reminder) const
    {
        if (this->layout_->centerings(field.quantity())[0] == this->layout_->dual)
        {
            if (reminder < 0.5)
                return field(iCell - 1) * (1.0 - reminder) + field(iCell) * reminder;
            else
                return field(iCell) * (1.0 - reminder) + field(iCell + 1) * reminder;
        }
        return field(iCell) * (1.0 - reminder) + field(iCell + 1) * reminder;
    }
};


#endif
