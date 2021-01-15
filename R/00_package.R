################################################################################################################################

# Name: freqpcr
# Short description: Interval estimation of population allele frequency based on quantitative PCR double-Delta Cq measures from bulk samples
# Author: Masaaki Sudo (NARO, JAPAN)
# Maintainer: Masaaki Sudo
## https://sudori.info/english.html
## https://orcid.org/0000-0001-9834-9857

# License:
# Copyright 2020 Masaaki Sudo
# The source code is distributed under the GNU General Public License version 3 (GPLv3)

## This program is free software: you can redistribute it and/or modify it under the terms of
## the GNU General Public License as published by the Free Software Foundation, either version 3 of the License,
## or (at your option) any later version.
## This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
## You should have received a copy of the GNU General Public License along with this program.
## If not, see <http://www.gnu.org/licenses/>.

#' The freqpcr package
#'
#' Maximum Interval estimation of population allele frequency based on quantitative PCR double-Delta Cq measures from bulk samples.
#'
#' @docType package
#' @importFrom methods new
#' @importFrom cubature cubintegrate
#' @importFrom stats dbeta dbinom dgamma dmultinom dnorm nlm optim plogis qlogis qnorm rbinom rgamma rmultinom rnorm
#' @importFrom utils flush.console
#' @keywords internal
#' @encoding UTF-8
#' @aliases freqpcr-package
"_PACKAGE"
