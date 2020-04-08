import math

R = 6378137
K0 = 0.9996
E = 0.00669438
E2 = math.pow(E, 2)
E3 = math.pow(E, 3)
E_P2 = E / (1 - E)
SQRT_E = math.sqrt(1 -E)
_E = (1 - SQRT_E) / (1 + SQRT_E)
_E2 = math.pow(_E, 2)
_E3 = math.pow(_E, 3)
_E4 = math.pow(_E, 4)
_E5 = math.pow(_E, 5)
M1 = 1 - E / 4 - 3 * E2 / 64 - 5 * E3 / 256
M2 = 3 * E / 8 + 3 * E2 / 32 + 45 * E3 / 1024
M3 = 15 * E2 / 256 + 45 * E3 / 1024
M4 = 35 * E3 / 3072
P2 = 3 / 2 * _E - 27 / 32 * _E3 + 269 / 512 * _E5
P3 = 21 / 16 * _E2 - 55 / 32 * _E4
P4 = 151 / 96 * _E3 - 417 / 128 * _E5
P5 = 1097 / 512 * _E4
MRHO = 206264.8062
PI2 = 6.283185307179586

ZONE_LETTERS = "CDEFGHJKLMNPQRSTUVWXX"

WGS84_ELLIPSOID = {
  "a": 6378137.0,
  "alpha": 1 / 298.257223563,
  "e2": 0.0066943799901413165
}

SK42_ELLIPSOID = {
  "a": 6378245.0,
  "alpha": 1 / 298.3,
  "e2": 0.006693421622965943
}

WGS84_SK42_TRANSFORM = {
  "dX": -23.92,
  "dY": 141.27,
  "dZ": 80.9,
  "omegaX": 0,
  "omegaY": -0.0000016968478842708164,
  "omegaZ": -0.000003975472186005913,
  "m": 0.12e-6
}

SK42_WGS84_TRANSFORM = {
  "dX": 23.92,
  "dY": -141.27,
  "dZ": -80.9,
  "omegaX": 0,
  "omegaY": 0.0000016968478842708164,
  "omegaZ": 0.000003975472186005913,
  "m": -0.12e-6
}

def to_deg(rad):
  return rad / math.pi * 180

def to_rad(deg):
  return deg * math.pi / 180

"""
get zone letter from latitude number
"""
def lat_to_zone_letter(latitude):
  if -80 <= latitude and latitude <= 84:
    return ZONE_LETTERS[math.floor((latitude + 80) / 8)]
  else:
    return None

"""
get zone number from lat/lon numbers
"""
def latlon_to_zone_number(latitude, longitude):
  if 56 <= latitude and latitude < 64 and 3 <= longitude and longitude < 12:
    return 32

  if 72 <= latitude and latitude <= 84 and longitude >= 0:
    if longitude < 9: return 31
    if longitude < 21: return 33
    if longitude < 33: return 35
    if longitude < 42: return 37

  return math.floor((longitude + 180) / 6) + 1


"""
convert zone number to central longitude
"""
def zone_number_to_central_lon(zn):
  return (zn - 1) * 6 - 180 + 3

"""
convert wgs84 props to utm object
"""
def latlon_utm(latitude, longitude):
  latRad = to_rad(latitude)
  latSin = math.sin(latRad)
  latCos = math.cos(latRad)

  latTan = math.tan(latRad)
  latTan2 = math.pow(latTan, 2)
  latTan4 = math.pow(latTan, 4)

  zoneNum = latlon_to_zone_number(latitude, longitude)
  zoneLetter = lat_to_zone_letter(latitude)

  lonRad = to_rad(longitude)
  centralLon = zone_number_to_central_lon(zoneNum)
  centralLonRad = to_rad(centralLon)

  n = R / math.sqrt(1 - E * latSin * latSin)
  c = E_P2 * latCos * latCos

  a = latCos * (lonRad - centralLonRad)
  a2 = math.pow(a, 2)
  a3 = math.pow(a, 3)
  a4 = math.pow(a, 4)
  a5 = math.pow(a, 5)
  a6 = math.pow(a, 6)

  m = R * (M1 * latRad - M2 * math.sin(2 * latRad) + M3 * math.sin(4 * latRad) - M4 * math.sin(6 * latRad))
  easting  = K0 * n * (a + a3 / 6 * (1 - latTan2 + c) + a5 / 120 * (5 - 18 * latTan2 + latTan4 + 72 * c - 58 * E_P2)) + 500000
  northing = K0 * (m + n * latTan * (a2 / 2 + a4 / 24 * (5 - latTan2 + 9 * c + 4 * c * c) + a6 / 720 * (61 - 58 * latTan2 + latTan4 + 600 * c - 330 * E_P2)))

  if latitude < 0:
    northing += 1e7

  return {
    "easting": easting,
    "northing": northing,
    "zoneNum": zoneNum,
    "zoneLetter": zoneLetter
  }

"""
convert utm props to wgs84 object
"""
def utm_latlon(easting, northing, zoneNum, zoneLetter, northern):
  if zoneLetter:
    zoneLetter = str.upper(zoneLetter)
    northern = zoneLetter >= "N"

  x = easting - 500000.0
  y = northing

  if not northern: y -= 1e7

  m = y / K0
  mu = m / (R * M1)

  pRad = mu + \
    P2 * math.sin(2 * mu) + \
    P3 * math.sin(4 * mu) + \
    P4 * math.sin(6 * mu) + \
    P5 * math.sin(8 * mu)

  pSin = math.sin(pRad)
  pSin2 = math.pow(pSin, 2)
  pCos = math.cos(pRad)
  pTan = math.tan(pRad)
  pTan2 = math.pow(pTan, 2)
  pTan4 = math.pow(pTan, 4)

  epSin = 1 - E * pSin2
  epSinSqrt = math.sqrt(epSin)

  n = R / epSinSqrt
  r = (1 - E) / epSin
  c = _E * pCos * pCos
  c2 = c * c

  d = x / (n * K0)
  d2 = math.pow(d, 2)
  d3 = math.pow(d, 3)
  d4 = math.pow(d, 4)
  d5 = math.pow(d, 5)
  d6 = math.pow(d, 6)

  latitude = pRad - (pTan / r) * (d2 / 2 - d4 / 24 * (5 + 3 * pTan2 + 10 * c - 4 * c2 - 9 * E_P2)) + d6 / 720 * (61 + 90 * pTan2 + 298 * c + 45 * pTan4 - 252 * E_P2 - 3 * c2)
  longitude = (d - d3 / 6 * (1 + 2 * pTan2 + c) + d5 / 120 * (5 - 2 * c + 28 * pTan2 - 3 * c2 + 8 * E_P2 + 24 * pTan4)) / pCos

  return {
    "latitude": to_deg(latitude),
    "longitude": to_deg(longitude) + zone_number_to_central_lon(zoneNum)
  }

"""
convert wgs84 props to mgrs string
"""
def latlon_mgrs(Lat, Long):
  if Lat < -80: return "Too far South"
  if Lat > 84: return "Too far North"

  c = 1 + math.floor((Long + 180) / 6)
  e = c * 6 - 183
  k = Lat * math.pi / 180
  l = Long * math.pi / 180
  m = e * math.pi / 180
  n = math.cos(k)
  o = 0.006739496819936062 * math.pow(n, 2)
  p = 40680631590769 / (6356752.314 * math.sqrt(1 + o))
  q = math.tan(k)
  r = q * q
  s = (r * r * r) - math.pow(q, 6)
  t = l - m
  u = 1.0 - r + o
  v = 5.0 - r + 9 * o + 4.0 * (o * o)
  w = 5.0 - 18.0 * r + (r * r) + 14.0 * o - 58.0 * r * o
  x = 61.0 - 58.0 * r + (r * r) + 270.0 * o - 330.0 * r * o
  y = 61.0 - 479.0 * r + 179.0 * (r * r) - (r * r * r)
  z = 1385.0 - 3111.0 * r + 543.0 * (r * r) - (r * r * r)
  aa = p * n * t + (p / 6.0 * math.pow(n, 3) * u * math.pow (t, 3)) + (p / 120.0 * math.pow(n, 5) * w * math.pow (t, 5)) + (p / 5040.0 * math.pow(n, 7) * y * math.pow(t, 7))
  ab = 6367449.14570093 * (k - (0.00251882794504 * math.sin(2 * k)) + (0.00000264354112 * math.sin(4 * k)) - (0.00000000345262 * math.sin(6 * k)) + (0.000000000004892 * math.sin(8 * k))) + (q / 2.0 * p * math.pow(n,2) * math.pow (t,2)) + (q / 24.0 * p * math.pow(n,4) * v * math.pow(t,4)) + (q / 720.0 * p * math.pow(n,6) * x * math.pow(t,6)) + (q / 40320.0 * p * math.pow(n,8) * z * math.pow(t,8))

  aa = aa * 0.9996 + 500000.0
  ab = ab * 0.9996

  if ab < 0.0: ab += 10000000.0

  ad = ZONE_LETTERS[math.floor(Lat / 8 + 10)]
  ae = math.floor(aa / 100000)
  af = ['ABCDEFGH', 'JKLMNPQR', 'STUVWXYZ'][(c - 1) % 3][ae - 1]
  ag = math.floor(ab / 100000) % 20
  ah = ['ABCDEFGHJKLMNPQRSTUV', 'FGHJKLMNPQRSTUVABCDE'][(c - 1) % 2][ag]

  def pad(val):
    if val < 10:
      val = "0000" + val
    elif val < 100:
      val = "000" + val
    elif val < 1000:
      val = "00" + val
    elif val < 10000:
      val = "0" + val

    return val

  aa = math.floor(aa % 100000)
  aa = pad(aa)
  ab = math.floor(ab % 100000)
  ab = pad(ab)

  return str(c) + str(ad) + " " + str(af) + str(ah) + " " + str(aa) + " " + str(ab)

"""
convert mgrs string to wgs84 object
"""
def mgrs_latlon(mgrs_string):
  b = mgrs_string.split()
  
  if len(b) != 4:
    return None

  c = b[0][0] if len(b[0]) < 3 else b[0][0:2]
  d = b[0][1] if len(b[0]) < 3 else b[0][2]

  e = (int(c) * 6 - 183) * math.pi / 180
  f = ["ABCDEFGH", "JKLMNPQR", "STUVWXYZ"][(int(c) - 1) % 3].index(b[1][0]) + 1
  g = ZONE_LETTERS.index(d)
  h = ["ABCDEFGHJKLMNPQRSTUV", "FGHJKLMNPQRSTUVABCDE"][(int(c) - 1) % 2].index(b[1][1])
  i = [1.1, 2.0, 2.8, 3.7, 4.6, 5.5, 6.4, 7.3, 8.2, 9.1, 0, 0.8, 1.7, 2.6, 3.5, 4.4, 5.3, 6.2, 7.0, 7.9]
  j = [0, 2, 2, 2, 4, 4, 6, 6, 8, 8, 0, 0, 0, 2, 2, 4, 4, 6, 6, 6]
  k = i[g]
  l = int(j[g]) + h / 10

  if l < k: l += 2

  m = f * 100000.0 + int(b[2])
  n = l * 1000000 + int(b[3])

  m -= 500000.0

  if d < "N": n -= 10000000.0

  m /= 0.9996
  n /= 0.9996

  o = n / 6367449.14570093
  p = o + (0.0025188266133249035 * math.sin(2.0 * o)) + (0.0000037009491206268 * math.sin(4.0 * o)) + (0.0000000074477705265 * math.sin(6.0 * o)) + (0.0000000000170359940 * math.sin(8.0 * o))

  q = math.tan(p)
  r = q * q
  s = r * r
  t = math.cos(p)
  u = 0.006739496819936062 * math.pow(t, 2)
  v = 40680631590769 / (6356752.314 * math.sqrt(1 + u))
  w = v
  x = 1.0 / (w * t); w *= v
  y = q / (2.0 * w); w *= v
  z = 1.0 / (6.0 * w * t); w *= v
  aa = q / (24.0 * w); w *= v
  ab = 1.0 / (120.0 * w * t); w *= v
  ac = q / (720.0 * w); w *= v
  ad = 1.0 / (5040.0 * w * t); w *= v
  ae = q / (40320.0 * w)
  af = -1.0 - u
  ag = -1.0 - 2 * r - u
  ah = 5.0 + 3.0 * r + 6.0 * u - 6.0 * r * u - 3.0 * (u * u) - 9.0 * r * (u * u)
  ai = 5.0 + 28.0 * r + 24.0 * s + 6.0 * u + 8.0 * r * u
  aj = -61.0 - 90.0 * r - 45.0 * s - 107.0 * u + 162.0 * r * u
  ak = -61.0 - 662.0 * r - 1320.0 * s - 720.0 * (s * r)
  al = 1385.0 + 3633.0 * r + 4095.0 * s + 1575 * (s * r)
  lat = p + y * af * (m * m) + aa * ah * math.pow(m,4) + ac * aj * math.pow(m,6) + ae * al * math.pow(m,8)
  lng = e + x * m + z * ag * math.pow(m, 3) + ab * ai * math.pow(m, 5) + ad * ak * math.pow(m, 7)

  lat = lat * 180 / math.pi
  lng = lng * 180 / math.pi

  return {
    "latitude": lat,
    "longitude": lng
  }


def elliptic_transform(latlon, ellips1, ellips2, tr):
  a = (ellips2["a"] + ellips1["a"]) / 2
  dA = ellips2["a"] - ellips1["a"]
  e2  = (ellips2["e2"] + ellips1["e2"]) / 2
  de2 = ellips2["e2"] - ellips1["e2"]

  B = latlon[0]
  L1 = latlon[1]
  H = 0

  M = a * (1 - e2) * math.pow(1 - e2 * math.pow(math.sin(B), 2), -1.5)
  N = a * math.pow(1 - e2 * math.pow(math.sin(B), 2), -0.5)

  dB = (MRHO / (M + H)) * ((N / a) * e2 * math.sin(B) * math.cos(B) * dA + (math.pow(N, 2) / math.pow(a, 2) + 1) * N * math.sin(B) * math.cos(B) * de2) / 2 - (tr.dX * math.cos(L1) + tr.dY * math.sin(L1) * math.sin(B) + tr.dZ * math.cos(B)) - tr.omegaX * math.sin(L1) * (1 + e2 * math.cos(2 * B)) + tr.omegaY * math.cos(L1) * (1 + e2 * math.cos(2 * B)) - MRHO * tr.m * e2 * math.sin(B) * math.cos(B)
  dL = (MRHO / ((N + H) * math.cos(B))) * (-tr.dX * math.sin(L1) + tr["dY"] * math.cos(L1)) + math.tan(B) *(1 - e2) * (tr["omegaX"] * math.cos(L1) + tr["omegaY"] * math.sin(L1)) - tr["omegaZ"]
  
  
  return [
    latlon[0] + dB / 3600,
    latlon[1] + dL / 3600
  ]

"""
convert longitude to pixel x
"""
def lon_x(lon, zoom):
  return math.floor((lon + 180) / 360 * math.pow(2, zoom))

"""
convert latitude to pixel y
"""
def lat_y(lat, zoom):
  return math.floor((1 - math.log(math.tan(lat * math.pi / 180) + 1 / math.cos(lat * math.pi / 180)) / math.pi) / 2 * math.pow(2, zoom))


                                    
                          
