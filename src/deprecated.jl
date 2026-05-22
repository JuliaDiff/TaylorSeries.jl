# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Luis Benet & David P. Sanders
# UNAM
#
# MIT Expat license
#

Base.@deprecate set_variables(args...; kwargs...) variables!(args...; kwargs...)
Base.@deprecate get_variables(args...; kwargs...) variables(args...; kwargs...)
