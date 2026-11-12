"""Minimal vector PDF plotting (no external dependencies).

Coordinate system: points (1/72 inch), origin at the lower-left corner.
Only the features needed by the paper figures are implemented: frames,
tick marks, labels, polygonal lines (solid / dashed / dotted), markers,
and a simple legend.
"""

import math, zlib


def _esc(s):
    return s.replace("\\", r"\\").replace("(", r"\(").replace(")", r"\)")


class Canvas:
    def __init__(self, width, height):
        self.w = width
        self.h = height
        self.ops = []
        self._data = None
        self._box = None
        self.images = []

    # ---------------- basic graphics ----------------
    def _out(self, s):
        self.ops.append(s)

    def axes(self, x0, y0, x1, y1):
        """reserve a plotting rectangle (figure points)"""
        self._box = (x0, y0, x1, y1)

    def data_range(self, xmin, xmax, ymin, ymax, logx=False, logy=False):
        self._data = (xmin, xmax, ymin, ymax, logx, logy)

    def _tf(self, x, y):
        xmin, xmax, ymin, ymax, logx, logy = self._data
        if logx:
            x = math.log10(x)
            xmin, xmax = math.log10(xmin), math.log10(xmax)
        if logy:
            y = math.log10(y)
            ymin, ymax = math.log10(ymin), math.log10(ymax)
        x0, y0, x1, y1 = self._box
        px = x0 + (x - xmin) / (xmax - xmin) * (x1 - x0)
        py = y0 + (y - ymin) / (ymax - ymin) * (y1 - y0)
        return px, py

    def frame(self, lw=0.8):
        x0, y0, x1, y1 = self._box
        self._out("%.2f w" % lw)
        self._out("%.2f %.2f %.2f %.2f re S" % (x0, y0, x1 - x0, y1 - y0))

    def polyline(self, pts, lw=0.8, dash=None, color=(0, 0, 0)):
        if len(pts) == 0:
            return
        self._out("%.3f %.3f %.3f RG" % color)
        self._out("%.2f w" % lw)
        if dash:
            self._out("[%s] 0 d" % " ".join("%.2f" % d for d in dash))
        else:
            self._out("[] 0 d")
        p0 = self._tf(*pts[0])
        self._out("%.2f %.2f m" % p0)
        for p in pts[1:]:
            q = self._tf(*p)
            self._out("%.2f %.2f l" % q)
        self._out("S")

    def marker(self, x, y, kind="circle", size=2.0, color=(0, 0, 0), filled=True):
        px, py = self._tf(x, y)
        self._out("%.3f %.3f %.3f RG" % color)
        self._out("%.3f %.3f %.3f rg" % color)
        self._out("0.8 w")
        if kind == "circle":
            k = 0.5523 * size
            self._out("%.2f %.2f m" % (px + size, py))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px + size, py + k, px + k, py + size, px, py + size))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px - k, py + size, px - size, py + k, px - size, py))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px - size, py - k, px - k, py - size, px, py - size))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px + k, py - size, px + size, py - k, px + size, py))
            self._out("h " + ("f" if filled else "S"))
        elif kind == "square":
            self._out("%.2f %.2f %.2f %.2f re %s" % (px - size, py - size, 2 * size, 2 * size,
                                                     "f" if filled else "S"))
        elif kind == "triangle":
            self._out("%.2f %.2f m" % (px, py + size * 1.2))
            self._out("%.2f %.2f l" % (px + size * 1.1, py - size * 0.9))
            self._out("%.2f %.2f l" % (px - size * 1.1, py - size * 0.9))
            self._out("h " + ("f" if filled else "S"))
        elif kind == "diamond":
            self._out("%.2f %.2f m" % (px, py + size * 1.2))
            self._out("%.2f %.2f l" % (px + size * 1.2, py))
            self._out("%.2f %.2f l" % (px, py - size * 1.2))
            self._out("%.2f %.2f l" % (px - size * 1.2, py))
            self._out("h " + ("f" if filled else "S"))

    def text(self, x, y, s, size=8.0, align="left", font="F1", angle=0.0):
        w = 0.5 * size * len(s) if font in ("F1", "F2") else 0.45 * size * len(s)
        if align == "center":
            x = x - w / 2
        elif align == "right":
            x = x - w
        self._out("BT /%s %.2f Tf" % (font, size))
        if angle:
            a = math.radians(angle)
            self._out("%.4f %.4f %.4f %.4f %.2f %.2f Tm" %
                      (math.cos(a), math.sin(a), -math.sin(a), math.cos(a), x, y))
        else:
            self._out("1 0 0 1 %.2f %.2f Tm" % (x, y))
        self._out("(%s) Tj ET" % _esc(s))

    # ---------------- graphics in figure coordinates ----------------
    def fig_line(self, x0, y0, x1, y1, lw=0.8, dash=None, color=(0, 0, 0)):
        self._out("%.3f %.3f %.3f RG" % color)
        self._out("%.2f w" % lw)
        if dash:
            self._out("[%s] 0 d" % " ".join("%.2f" % d for d in dash))
        else:
            self._out("[] 0 d")
        self._out("%.2f %.2f m %.2f %.2f l S" % (x0, y0, x1, y1))

    def fig_marker(self, px, py, kind="circle", size=2.0, filled=True):
        self._out("0 0 0 RG")
        self._out("0 0 0 rg")
        self._out("0.8 w")
        if kind == "circle":
            k = 0.5523 * size
            self._out("%.2f %.2f m" % (px + size, py))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px + size, py + k, px + k, py + size, px, py + size))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px - k, py + size, px - size, py + k, px - size, py))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px - size, py - k, px - k, py - size, px, py - size))
            self._out("%.2f %.2f %.2f %.2f %.2f %.2f c" % (px + k, py - size, px + size, py - k, px + size, py))
            self._out("h " + ("f" if filled else "S"))
        elif kind == "square":
            self._out("%.2f %.2f %.2f %.2f re %s" % (px - size, py - size, 2 * size, 2 * size,
                                                     "f" if filled else "S"))
        elif kind == "triangle":
            self._out("%.2f %.2f m" % (px, py + size * 1.2))
            self._out("%.2f %.2f l" % (px + size * 1.1, py - size * 0.9))
            self._out("%.2f %.2f l" % (px - size * 1.1, py - size * 0.9))
            self._out("h " + ("f" if filled else "S"))
        elif kind == "diamond":
            self._out("%.2f %.2f m" % (px, py + size * 1.2))
            self._out("%.2f %.2f l" % (px + size * 1.2, py))
            self._out("%.2f %.2f l" % (px, py - size * 1.2))
            self._out("%.2f %.2f l" % (px - size * 1.2, py))
            self._out("h " + ("f" if filled else "S"))

    def sci_label(self, x, y, expo, size=7.0):
        """tick label  10^expo  (x,y = right end of the label)"""
        w1 = 0.5 * size * 2
        we = 0.5 * size * 0.68 * len(expo)
        x = x - (w1 + we + 0.5)
        self.text(x, y, "10", size)
        self.text(x + w1 + 0.4, y + 0.40 * size, expo, size * 0.68)

    def ticks(self, xticks, yticks, xlabels, ylabels, size=7.5, length=2.5,
              xlabel=None, ylabel=None, minor=0):
        x0, y0, x1, y1 = self._box
        for xv, lab in zip(xticks, xlabels):
            px, py = self._tf(xv, y0)
            self._out("0.6 w")
            self._out("%.2f %.2f m %.2f %.2f l S" % (px, y0, px, y0 - length))
            self.text(px, y0 - length - size - 1.5, lab, size, "center")
        for yv, lab in zip(yticks, ylabels):
            px, py = self._tf(x0, yv)
            self._out("0.6 w")
            self._out("%.2f %.2f m %.2f %.2f l S" % (x0, py, x0 - length, py))
            self.text(x0 - length - 2.5, py - 0.35 * size, lab, size, "right")
        if xlabel:
            self.text((x0 + x1) / 2, y0 - 3.2 * size - 4, xlabel, size + 1, "center")
        if ylabel:
            self.text(x0 - 4.2 * size, (y0 + y1) / 2, ylabel, size + 1, "center", angle=90)

    def image(self, x, y, w, h, rgb, iw, ih):
        """draw an RGB image (bytes, iw x ih) into the rectangle (x,y,w,h)"""
        self.images.append((rgb, iw, ih))
        k = len(self.images)
        self._out("q %.2f 0 0 %.2f %.2f %.2f cm /Im%d Do Q" % (w, h, x, y, k))

    def save(self, path):
        content = "\n".join(self.ops) + "\n"
        objs = []
        objs.append("<< /Type /Catalog /Pages 2 0 R >>")
        objs.append("<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
        xobj = ""
        if self.images:
            xobj = " /XObject << " + " ".join("/Im%d %d 0 R" % (k, 7 + k)
                                              for k in range(1, len(self.images) + 1)) + " >>"
        objs.append("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] "
                    "/Resources << /Font << /F1 5 0 R /F2 6 0 R /F3 7 0 R >>%s >> "
                    "/Contents 4 0 R >>" % (self.w, self.h, xobj))
        objs.append("<< /Length %d >>\nstream\n%s\nendstream" % (len(content), content))
        objs.append("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
        objs.append("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>")
        objs.append("<< /Type /Font /Subtype /Type1 /BaseFont /Times-Roman >>")
        for rgb, iw, ih in self.images:
            comp = zlib.compress(rgb, 9)
            objs.append("<< /Type /XObject /Subtype /Image /Width %d /Height %d "
                        "/ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode "
                        "/Length %d >>\nstream\n" % (iw, ih, len(comp))
                        + comp.decode("latin-1") + "\nendstream")

        out = b"%PDF-1.4\n"
        offsets = []
        for i, o in enumerate(objs, start=1):
            offsets.append(len(out))
            out += ("%d 0 obj\n%s\nendobj\n" % (i, o)).encode("latin-1")
        xref = len(out)
        out += ("xref\n0 %d\n" % (len(objs) + 1)).encode()
        out += b"0000000000 65535 f \n"
        for off in offsets:
            out += ("%010d 00000 n \n" % off).encode()
        out += ("trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n"
                % (len(objs) + 1, xref)).encode()
        open(path, "wb").write(out)
