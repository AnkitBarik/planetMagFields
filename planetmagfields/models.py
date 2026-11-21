#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Registry of planets and field models.

Each planet entry holds its default model and the file format of its data files
("generic", "igrf", "jrm" or "axisymmetric"). Model entries can override the
format and lmax and carry citation information.

Coordinate convention: for every planet the Gauss coefficients are for Schmidt
semi-normalised Legendre functions without the Condon-Shortley phase, and the
longitude phi is east longitude in the body-fixed frame of the field model
(e.g. geographic longitude for Earth, System III for Jupiter).
"""

PLANETS = {
    "mercury":  {"default_model": "wardinski2019"},
    "earth":    {"default_model": "igrf14", "format": "igrf"},
    "jupiter":  {"default_model": "jrm33"},
    "saturn":   {"default_model": "cassini11+", "format": "axisymmetric"},
    "uranus":   {"default_model": "herbert2009"},
    "neptune":  {"default_model": "holme1996"},
    "ganymede": {"default_model": "kivelson2002"},
}

MODELS = {
    "anderson2012":  {"name": "Anderson et al. 2012",
                      "url": "https://doi.org/10.1029/2012JE004159",
                      "format": "axisymmetric"},
    "thebault2018":  {"name": "Thébault et al. 2018",
                      "url": "https://doi.org/10.1016/j.pepi.2017.07.001"},
    "wardinski2019": {"name": "Wardinski et al. 2019",
                      "url": "https://doi.org/10.1029/2018JE005835"},
    "igrf14":        {"name": "IGRF 14",
                      "url": "https://doi.org/10.1186/s40623-020-01288-x",
                      "data_url": "https://doi.org/10.5281/zenodo.14012303"},
    "vip4":          {"name": "VIP4 (Connerney et al. 1998)",
                      "url": "https://doi.org/10.1029/97JA03726"},
    "jrm09":         {"name": "JRM09 (Connerney et al. 2018)",
                      "url": "https://doi.org/10.1002/2018GL077312",
                      "format": "jrm"},
    "jrm33":         {"name": "JRM33 (Connerney et al. 2022)",
                      "url": "https://doi.org/10.1029/2021JE007055",
                      "format": "jrm", "lmax": 18},
    "cassinisoi":    {"name": "Cassini SOI (Burton et al. 2009)",
                      "url": "https://doi.org/10.1016/j.pss.2009.04.008"},
    "cassini11":     {"name": "Cassini 11 (Dougherty et al. 2018)",
                      "url": "https://doi.org/10.1126/science.aat5434"},
    "cassini11+":    {"name": "Cassini 11+ (Cao et al. 2020)",
                      "url": "https://doi.org/10.1016/j.icarus.2019.113541"},
    "connerney1987": {"name": "Connerney et al. 1987",
                      "url": "https://doi.org/10.1029/JA092iA13p15329"},
    "holme1996":     {"name": "Holme & Bloxham 1996",
                      "url": "https://doi.org/10.1029/95JE03437"},
    "herbert2009":   {"name": "Herbert 2009",
                      "url": "https://doi.org/10.1029/2009JA014394"},
    "connerney1991": {"name": "Connerney et al. 1991",
                      "url": "https://doi.org/10.1029/91JA01165"},
    "kivelson2002":  {"name": "Kivelson et al. 2002",
                      "url": "https://doi.org/10.1006/icar.2002.6834"},
}

planetlist = list(PLANETS)


def check_planet(planetname):
    planetname = planetname.lower()
    if planetname not in PLANETS:
        raise ValueError("Unknown planet '%s', must be one of %s"
                         % (planetname, planetlist))
    return planetname


def default_model(planetname):
    return PLANETS[check_planet(planetname)]["default_model"]


def model_info(model):
    return MODELS.get(model.lower(), {})


def model_format(planetname, model):
    fmt = model_info(model).get("format")
    return fmt or PLANETS[check_planet(planetname)].get("format", "generic")


def phase_factor(planetname, m):
    """Factor that removes the Condon-Shortley phase included in scipy's
    spherical harmonics. Gauss coefficients of all planets are defined
    without it, so this is (-1)**m for every planet; planetname is kept for
    backward compatibility."""
    return (-1)**m
