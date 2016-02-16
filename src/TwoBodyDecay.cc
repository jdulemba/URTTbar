#include "TwoBodyDecay.h"
#include "Permutation.h"
#include "GenObject.h"

TwoBodyDecay::TwoBodyDecay(Permutation &p):
  TwoBodyDecay(
    new TwoBodyDecay(
      new TwoBodyDecay(p.BLep()),
      new TwoBodyDecay(p.L(), p.NuPtr())
      ),
    new TwoBodyDecay(
      new TwoBodyDecay(p.BHad()),
      new TwoBodyDecay(p.WJa(), p.WJb())
      )
    ) {}

TwoBodyDecay::TwoBodyDecay(GenTTBar &g):
  TwoBodyDecay(
    new TwoBodyDecay(
      &g.top,
      new TwoBodyDecay(g.top.b),
      new TwoBodyDecay(g.top.W.first, g.top.W.second)
      ),
    new TwoBodyDecay(
      &g.tbar,
      new TwoBodyDecay(g.tbar.b),
      new TwoBodyDecay(g.tbar.W.first, g.top.W.second)
      )
    )  {}       

double TwoBodyDecay::cosThetaStar_decay() {
  TVector3 direction = unit3D();
  TwoBodyDecay cmframe = to_CM();
  TVector3 fst_dir = cmframe.fst()->unit3D();
  return direction.Dot(fst_dir);
}
// 
// 
// 
// 
// 
// 
// 
// 
// 
