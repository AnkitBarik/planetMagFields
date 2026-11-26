#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Registry of planets and field models, loaded from data/models.json.

The file holds one entry per planet under "planets"; their order sets the order
of ``planetlist``. Each planet entry has

  "default_model" : model used when none is given
  "format"        : layout of its data files ("generic", "igrf", "jrm" or
                    "axisymmetric"); "generic" if omitted
  "models"        : one entry per available model, with citation information
                    ("name", "url" and optionally "data_url"). A model entry
                    can override the planet's "format" and the "lmax" read
                    from the data file.

To add a model, place its data file in data/ as <planet>_<model>.dat and add
an entry for <model> under that planet's "models" in data/models.json.

Coordinate convention: for every planet the Gauss coefficients are for Schmidt
semi-normalised Legendre functions without the Condon-Shortley phase, and the
longitude phi is east longitude in the body-fixed frame of the field model
(e.g. geographic longitude for Earth, System III for Jupiter).
"""

import json
import os
import sys

from .termcolors import red

_REGISTRY_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'data', 'models.json')

with open(_REGISTRY_FILE, encoding='utf-8') as _f:
    _registry = json.load(_f)

PLANETS = _registry["planets"]

# Flat {model: info} view of all planets' models, kept for backward
# compatibility. Model names shared between planets (e.g. holme1996 for Uranus
# and Neptune) map to the last planet's entry; use model_info(model, planet)
# for planet-specific information.
MODELS = {model: info for planet in PLANETS.values()
          for model, info in planet["models"].items()}

planetlist = list(PLANETS)


def check_planet(planetname):
    planetname = planetname.lower()
    if planetname not in PLANETS:
        raise ValueError("Unknown planet '%s'.\n%s"
                         % (planetname, red("Planets available: %s"
                                            % ", ".join(planetlist), sys.stderr)))
    return planetname


def default_model(planetname):
    return PLANETS[check_planet(planetname)]["default_model"]


def model_info(model, planetname=None):
    """Registry entry of a model (citation, format/lmax overrides), or {} if
    the model is not listed. Give planetname to look the model up under that
    planet; without it the flat MODELS view is used."""
    model = model.lower()
    if planetname is None:
        return MODELS.get(model, {})
    return PLANETS[check_planet(planetname)]["models"].get(model, {})


def model_format(planetname, model):
    fmt = model_info(model, planetname).get("format")
    return fmt or PLANETS[check_planet(planetname)].get("format", "generic")


def phase_factor(planetname, m):
    """Factor that removes the Condon-Shortley phase included in scipy's
    spherical harmonics. Gauss coefficients of all planets are defined
    without it, so this is (-1)**m for every planet; planetname is kept for
    backward compatibility."""
    return (-1)**m
