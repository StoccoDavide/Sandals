/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef TESTS_GERDTS_PROBLEMS_AD_HH
#define TESTS_GERDTS_PROBLEMS_AD_HH

// The DAE optimal control examples of M. Gerdts, Optimal Control of ODEs and DAEs, 2nd edition, De
// Gruyter, 2023, written once as templates (derivatives by CppAD), each in the formulation of the
// book and, where meaningful, index-reduced with the lost constraints as invariants:
//  - Gerdts/Academic.hh:   3.1.13 cost on the algebraic state, 7.1.13 bang-singular-bang, 7.1.6
//                          degenerate problem, 1.1.17/5.7.4 Kostyukova (unsmoothed and smoothed);
//  - Gerdts/Mechanical.hh: 1.1.19/4.6.5/6.1.9 GGL pendulum, 5.7.3 ski driver, 5.7.6 half-car chassis;
//  - Gerdts/Process.hh:    8.2.8 semi-batch chemical reactor;
//  - Gerdts/Flight.hh:     1.1.10/6.1.10 emergency landing, free and on a prescribed path;
//  - Gerdts/Stokes.hh:     1.0.4/1.1.11/3.1.12 distributed control of the discretized Stokes flow.
// The minimum energy (1.0.1, 3.1.11, 3.2.13) and pendulum (1.1.9) examples are in DAEProblemsAD.hh.

#include "Gerdts/Academic.hh"
#include "Gerdts/Mechanical.hh"
#include "Gerdts/Process.hh"
#include "Gerdts/Flight.hh"
#include "Gerdts/Stokes.hh"

#endif // TESTS_GERDTS_PROBLEMS_AD_HH
