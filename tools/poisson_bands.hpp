/* poisson_bands.hpp
 *
 * Copyright (c) 2019 Luigi Pertoldi
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <utility>
#include <cmath>
#include <stdexcept>

#include "TMath.h"
#include "TColor.h"
#include "TBox.h"
#include "TH1.h"

namespace poiband {

    /*
     * smallest_poisson_interval(prob_coverage, poisson_mean)
     *
     * Computes the smallest interval boundaries covering `prob_coverage` area of a
     * discrete Poisson distribution with mean `poisson_mean` Returns a
     * `std::pair<float,float>` holding the lower and upper range
     */
    std::pair<float, float> smallest_poisson_interval(double cov, double mu) {

        // sanity check
        if (cov > 1 or cov < 0 or mu < 0) throw std::runtime_error("smallest_poisson_interval(): bad input");

        // initialize lower and upper edges to something
        std::pair<float, float> res = {mu, mu};

        if (mu > 50) { // gaussian approximation os OK
            res = {
                std::round(mu + TMath::NormQuantile((1-cov)/2)*sqrt(mu))-0.5,
                std::round(mu - TMath::NormQuantile((1-cov)/2)*sqrt(mu))+0.5
            };
        }
        else { // do the computation
            // start from the mode, which is the integer part of the mean
            int mode = std::floor(mu);
            int l = mode, u = mode; // let's start from here
            double prob = TMath::PoissonI(mode, mu); // probability covered by the interval

            // check if we're undercovering
            while (prob < cov) {
                // compute probabilities of points just ouside interval
                double prob_u = TMath::PoissonI(u+1, mu);
                double prob_l = TMath::PoissonI(l > 0 ? l-1 : 0, mu);

                // we expand on the right if:
                //  - the lower edge is already at zero
                //  - the prob of the right point is higher than the left
                if (l == 0 or prob_u > prob_l) {
                    u++; // expand interval
                    prob += prob_u; // update coverage
                }
                // otherwhise we expand on the left
                else if (prob_u < prob_l) {
                    l--;
                    prob += prob_l;
                }
                // if prob_u == prob_l we expand on both sides
                else {
                    u++; l--;
                    prob += prob_u + prob_l;
                }
            }
            res = {l == 0 ? 0 : l-0.5, u+0.5};
        }
        return res;
    }

    bool col_defined = false;
    /*
     * draw_poisson_bands(poisson_mean, box_x_left_location, box_size, normalize_to_mean)
     *
     * draws TBoxes corresponding to 68, 95 and 98 coverage intervals for poisson
     * (discrete) distribution with mean `poisson_mean`. The location and size of
     * the box on the x-axis must be provided with the second and third arguments.
     * The boolean argument `residuals` can be used to normalize boxes to the value
     * of the `poisson_mean`. The histogram onto which the bands will be drawn can
     * be provided as a last argument, and the boxes will be clipped to its frame
     * size.
     */
    void draw_poisson_bands(double mu, double x_low, double x_size, bool residuals = false, double scale = 1., TH1* h = nullptr) {

        if (h != nullptr and h->GetDimension() != 1) {
            throw std::runtime_error("draw_poisson_bands(): only 1D histograms are supported");
        }

        int col_idx = 9000;
        if (!col_defined) {
            new TColor(col_idx  , 238./255, 136./255, 102./255, "tol-lig-orange");
            new TColor(col_idx+1, 238./255, 221./255, 136./255, "tol-lig-lightyellow");
            new TColor(col_idx+2, 187./255, 204./255,  51./255, "tol-lig-pear");
            col_defined = true;
        }

        // calculate smallest intervals
        auto sig1 = smallest_poisson_interval(0.682, mu);
        auto sig2 = smallest_poisson_interval(0.954, mu);
        auto sig3 = smallest_poisson_interval(0.997, mu);

        if (residuals) {
            if (mu != 0) {
                sig1.first /= mu; sig1.second /= mu;
                sig2.first /= mu; sig2.second /= mu;
                sig3.first /= mu; sig3.second /= mu;
            }
            else {
                sig1.first = sig1.second = 1;
                sig2.first = sig2.second = 1;
                sig3.first = sig3.second = 1;
            }
        }
        else {
            sig1.first *= scale; sig1.second *= scale;
            sig2.first *= scale; sig2.second *= scale;
            sig3.first *= scale; sig3.second *= scale;
        }

        // compute interval centers (for plotting)
        auto cent_b1 = (sig1.second + sig1.first)/2;
        auto cent_b2 = (sig2.second + sig2.first)/2;
        auto cent_b3 = (sig3.second + sig3.first)/2;

        auto xdw = x_low;
        auto xup = x_low + x_size;

        // do now draw bands outside histogram frame
        if (h != nullptr) {
            auto xc1 = gPad->GetUxmin(); auto xc2 = gPad->GetUxmax();
            auto yc1 = gPad->GetUymin(); auto yc2 = gPad->GetUymax();

            if (sig1.first  < yc1) sig1.first  = yc1;
            if (sig2.first  < yc1) sig2.first  = yc1;
            if (sig3.first  < yc1) sig3.first  = yc1;
            if (sig1.second > yc2) sig1.second = yc2;
            if (sig2.second > yc2) sig2.second = yc2;
            if (sig3.second > yc2) sig3.second = yc2;
            if (xdw < xc1) xdw = xc1;
            if (xup > xc2) xup = xc2;
        }

        auto box_b1 = new TBox(xdw, sig1.first, xup, sig1.second);
        auto box_b2 = new TBox(xdw, sig2.first, xup, sig2.second);
        auto box_b3 = new TBox(xdw, sig3.first, xup, sig3.second);

        box_b3->SetFillColor(col_idx);
        box_b2->SetFillColor(col_idx+1);
        box_b1->SetFillColor(col_idx+2);

        if (box_b3->GetY1() != box_b3->GetY2()) box_b3->Draw();
        if (box_b2->GetY1() != box_b2->GetY2()) box_b2->Draw();
        if (box_b1->GetY1() != box_b1->GetY2()) box_b1->Draw();

        return;
    }
}

// vim: ft=cpp
