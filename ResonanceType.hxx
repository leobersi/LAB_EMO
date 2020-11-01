#ifndef RESONANCE_TYPE_H
#define RESONANCE_TYPE_H
#include "ParticleType.hxx"
class ResonanceType : public ParticleType {
public:
ResonanceType(const char* name, double mass, int charge, double width);
double GetWidth() const override;
void Print() const override;
private:
const double fWidth_; //lenght of risonance
};
#endif