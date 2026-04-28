from __future__ import annotations

from .viz import Animate3D


class Scene3D(Animate3D):
    """
    General astronomy-facing alias for the main 3D scene builder.

    `Scene3D` is intentionally additive: it reuses the existing `Animate3D`
    implementation so current notebooks and scripts remain valid.
    """

    pass
