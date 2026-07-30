"""Astronomy-facing alias for the main Oviz scene builder."""

from __future__ import annotations

from .viz import Animate3D


class Scene3D(Animate3D):
    """Astronomy-facing alias for :class:`oviz.Animate3D`.

    ``Scene3D`` reuses the same implementation so existing ``Animate3D``
    notebooks and scripts remain valid.
    """

    pass
