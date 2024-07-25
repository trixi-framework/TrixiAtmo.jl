"""
    Physical constants
"""
const EARTH_RADIUS = 6.37122f6  # m
const EARTH_GRAVITATIONAL_ACCELERATION = 9.80616f0  # m/s²
const EARTH_ROTATION_RATE = 7.292f-5  # rad/s
const SECONDS_PER_DAY = 8.64f4

abstract type AbstractCovariantEquations2D{NVARS} <: Trixi.AbstractEquations{2, NVARS} end

export CovariantLinearAdvectionEquation2D
include("covariant_advection.jl")
