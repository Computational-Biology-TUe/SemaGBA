"""
    <SemaGBA: A system dynamics model of the Semaglutide-responsive Gut-Brain Axis. 
    A model of how the brain and semaglutide regulate appetite and weight.>
    Copyright (C) 2026  Vivan Kennis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

using Pkg
Pkg.activate(".")
Pkg.instantiate()

function main(experiments)
    if isempty(experiments)
        experiments = [1,2,3,4]
    else
        experiments = parse.(Int, experiments)
    end

    if 1 in experiments
        # run exp 1
    end

end


main(ARGS)

