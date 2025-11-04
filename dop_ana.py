"""
LEO-PNT DOP Analyzer
====================

A tool for analyzing Dilution of Precision (DOP) metrics for LEO-based
Position, Navigation, and Timing (PNT) constellations.
"""

import numpy as np
from datetime import datetime, timedelta
from sgp4.api import Satrec, jday
import warnings

warnings.filterwarnings('ignore')


class DOPAnalyzer:
    """
    Analyzer for computing Dilution of Precision metrics for satellite constellations.
    
    The DOP metrics quantify the geometric strength of satellite configurations:
    - HDOP: Horizontal Dilution of Precision (2D position accuracy)
    - VDOP: Vertical Dilution of Precision (altitude accuracy)  
    - PDOP: Position Dilution of Precision (3D position accuracy)
    - GDOP: Geometric Dilution of Precision (position + time accuracy)
    - TDOP: Time Dilution of Precision (time accuracy)
    
    Lower DOP values indicate better geometric configuration and accuracy.
    """
    
    # WGS84 Earth model constants
    EARTH_A = 6378137.0  # Semi-major axis (m)
    EARTH_E2 = 0.00669437999014  # First eccentricity squared
    
    def __init__(self, tle_data, mask_angles=[5, 20, 30]):
        """
        Initialize the DOP analyzer.
        
        Args:
            tle_data (str): TLE data string containing satellite orbital elements
            mask_angles (list): List of elevation mask angles in degrees
        """
        self.mask_angles = mask_angles
        self.satellites = self._parse_tle_data(tle_data)
        
    def _parse_tle_data(self, tle_data):
        """
        Parse TLE (Two-Line Element) data into satellite objects.
        
        TLE format consists of three lines per satellite:
        - Line 0: Satellite name
        - Line 1: Orbital elements (epoch, inclination, etc.)
        - Line 2: Additional orbital elements
        
        Args:
            tle_data (str): Multi-line string containing TLE data
            
        Returns:
            list: List of dictionaries with 'name' and 'satrec' keys
        """
        lines = tle_data.strip().split('\n')
        satellites = []
        
        # Clean and filter empty lines
        clean_lines = [line.strip() for line in lines if line.strip()]
        
        # Parse in groups of 3 lines
        for i in range(0, len(clean_lines), 3):
            if i + 2 < len(clean_lines):
                name = clean_lines[i]
                line1 = clean_lines[i + 1]
                line2 = clean_lines[i + 2]
                
                # Validate TLE format
                if not (line1.startswith('1 ') and line2.startswith('2 ')):
                    print(f"Warning: Invalid TLE format for satellite {name}")
                    continue
                
                try:
                    satellite = Satrec.twoline2rv(line1, line2)
                    satellites.append({
                        'name': name,
                        'satrec': satellite
                    })
                except Exception as e:
                    print(f"Error parsing satellite {name}: {e}")
                    
        return satellites
    
    def _compute_satellite_position(self, satellite, jd, fr):
        """
        Compute satellite position in ECEF coordinates using SGP4.
        
        Args:
            satellite (dict): Satellite dictionary with 'satrec' key
            jd (float): Julian date
            fr (float): Fractional day
            
        Returns:
            np.array or None: Position vector [x, y, z] in km, or None if error
        """
        error, r, v = satellite['satrec'].sgp4(jd, fr)
        if error == 0:
            return np.array(r)
        return None
    
    def _ecef_to_lla(self, ecef_pos):
        """
        Convert ECEF (Earth-Centered Earth-Fixed) to LLA (Lat, Lon, Alt).
        
        Uses iterative algorithm for geodetic latitude computation.
        
        Args:
            ecef_pos (np.array): ECEF position [x, y, z] in km
            
        Returns:
            tuple: (latitude, longitude, altitude) in degrees, degrees, km
        """
        x, y, z = ecef_pos * 1000  # Convert km to m
        
        # Longitude
        lon = np.arctan2(y, x)
        
        # Latitude (iterative solution for geodetic latitude)
        p = np.sqrt(x**2 + y**2)
        lat = np.arctan2(z, p * (1 - self.EARTH_E2))
        
        for _ in range(5):  # Usually converges in 2-3 iterations
            N = self.EARTH_A / np.sqrt(1 - self.EARTH_E2 * np.sin(lat)**2)
            h = p / np.cos(lat) - N
            lat = np.arctan2(z, p * (1 - self.EARTH_E2 * N / (N + h)))
        
        # Altitude
        N = self.EARTH_A / np.sqrt(1 - self.EARTH_E2 * np.sin(lat)**2)
        h = p / np.cos(lat) - N
        
        return np.degrees(lat), np.degrees(lon), h / 1000
    
    def _compute_elevation_azimuth(self, sat_pos, observer_pos):
        """
        Compute elevation and azimuth angles from observer to satellite.
        
        This transforms the satellite position from ECEF to local
        East-North-Up (ENU) coordinates relative to the observer.
        
        Args:
            sat_pos (np.array): Satellite ECEF position [x, y, z] in km
            observer_pos (tuple): Observer (lat, lon, alt) in degrees, degrees, km
            
        Returns:
            tuple: (elevation, azimuth) in degrees
        """
        lat, lon, alt = observer_pos
        lat_rad, lon_rad = np.radians(lat), np.radians(lon)
        
        # Observer position in ECEF
        N = self.EARTH_A / np.sqrt(1 - self.EARTH_E2 * np.sin(lat_rad)**2)
        obs_ecef = np.array([
            (N + alt * 1000) * np.cos(lat_rad) * np.cos(lon_rad),
            (N + alt * 1000) * np.cos(lat_rad) * np.sin(lon_rad),
            (N * (1 - self.EARTH_E2) + alt * 1000) * np.sin(lat_rad)
        ]) / 1000  # Convert to km
        
        # Vector from observer to satellite
        sat_vector = sat_pos - obs_ecef
        
        # ECEF to ENU transformation matrix
        sin_lat, cos_lat = np.sin(lat_rad), np.cos(lat_rad)
        sin_lon, cos_lon = np.sin(lon_rad), np.cos(lon_rad)
        
        T = np.array([
            [-sin_lon, cos_lon, 0],
            [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
            [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat]
        ])
        
        enu = T @ sat_vector
        east, north, up = enu
        
        # Compute elevation and azimuth
        elevation = np.degrees(np.arctan2(up, np.sqrt(east**2 + north**2)))
        azimuth = np.degrees(np.arctan2(east, north))
        
        if azimuth < 0:
            azimuth += 360
            
        return elevation, azimuth
    
    def compute_dop(self, observer_pos, timestamp, mask_angle):
        """
        Compute all DOP metrics for a given observer position and time.

        The geometry matrix H relates satellite line-of-sight vectors to
        position/time errors. DOP values are derived from (H^T H)^-1.

        Args:
            observer_pos (tuple): Observer (lat, lon, alt) in degrees, degrees, km
            timestamp (datetime): Time of computation
            mask_angle (float): Minimum elevation angle in degrees

        Returns:
            dict: Dictionary containing GDOP, PDOP, HDOP, VDOP, TDOP, and num_satellites
        """
        jd, fr = jday(timestamp.year, timestamp.month, timestamp.day,
                      timestamp.hour, timestamp.minute, timestamp.second)

        visible_satellites = []

        # Find all visible satellites
        for satellite in self.satellites:
            sat_pos = self._compute_satellite_position(satellite, jd, fr)
            if sat_pos is not None:
                elevation, azimuth = self._compute_elevation_azimuth(sat_pos, observer_pos)
                if elevation >= mask_angle:
                    visible_satellites.append({
                        'elevation': elevation,
                        'azimuth': azimuth,
                        'position': sat_pos
                    })

        # Need at least 4 satellites for 3D+time solution
        if len(visible_satellites) < 4:
            return {
                'GDOP': float('nan'),
                'PDOP': float('nan'),
                'HDOP': float('nan'),
                'VDOP': float('nan'),
                'TDOP': float('nan'),
                'num_satellites': len(visible_satellites)
            }

        # Build geometry matrix H: each row is [east, north, up, clock]
        H_rows = []
        for sat in visible_satellites:
            el_rad = np.radians(sat['elevation'])
            az_rad = np.radians(sat['azimuth'])

            dx = np.sin(az_rad) * np.cos(el_rad)  # East
            dy = np.cos(az_rad) * np.cos(el_rad)  # North
            dz = np.sin(el_rad)                    # Up
            dt = 1.0                               # Clock term

            H_rows.append([dx, dy, dz, dt])

        H = np.array(H_rows, dtype=float)

        # 1) Rank check (need rank 4 for 3D + clock)
        if np.linalg.matrix_rank(H) < 4:
            return {
                'GDOP': float('nan'), 'PDOP': float('nan'), 'HDOP': float('nan'),
                'VDOP': float('nan'), 'TDOP': float('nan'),
                'num_satellites': len(visible_satellites)
            }

        # 2) SVD-based conditioning on H (more robust than checking H^T H only)
        s = np.linalg.svd(H, compute_uv=False)
        if (not np.all(np.isfinite(s))) or (s.min() <= 0) or (s.min() / s.max() < 1e-6):
            return {
                'GDOP': float('nan'), 'PDOP': float('nan'), 'HDOP': float('nan'),
                'VDOP': float('nan'), 'TDOP': float('nan'),
                'num_satellites': len(visible_satellites)
            }

        try:
            HTH = H.T @ H
            Q = np.linalg.pinv(HTH)

            qxx, qyy, qzz, qtt = Q[0, 0], Q[1, 1], Q[2, 2], Q[3, 3]
            HDOP = float(np.sqrt(qxx + qyy))
            VDOP = float(np.sqrt(qzz))
            PDOP = float(np.sqrt(qxx + qyy + qzz))
            TDOP = float(np.sqrt(qtt))
            GDOP = float(np.sqrt(qxx + qyy + qzz + qtt))

            # 3) Sanity gate: treat absurd values as invalid
            if (not np.isfinite(GDOP)) or (GDOP > 50.0) or (PDOP > 50.0):
                return {
                    'GDOP': float('nan'), 'PDOP': float('nan'), 'HDOP': float('nan'),
                    'VDOP': float('nan'), 'TDOP': float('nan'),
                    'num_satellites': len(visible_satellites)
                }

            return {
                'GDOP': GDOP, 'PDOP': PDOP, 'HDOP': HDOP,
                'VDOP': VDOP, 'TDOP': TDOP, 'num_satellites': len(visible_satellites)
            }

        except np.linalg.LinAlgError:
            return {
                'GDOP': float('nan'), 'PDOP': float('nan'), 'HDOP': float('nan'),
                'VDOP': float('nan'), 'TDOP': float('nan'),
                'num_satellites': len(visible_satellites)
            }
    
    def generate_global_dop_grid(self, timestamp, mask_angle, resolution=10):
        """
        Generate a global grid of DOP values.
        
        Args:
            timestamp (datetime): Time of computation
            mask_angle (float): Minimum elevation angle in degrees
            resolution (int): Grid spacing in degrees (smaller = more points)
            
        Returns:
            tuple: (lats, lons, dop_values) where dop_values is list of dicts
        """
        lats = []
        lons = []
        dop_values = []
        
        total_points = ((180 // resolution) + 1) * ((360 // resolution) + 1)
        print(f"Computing DOP for ~{total_points} grid points (resolution={resolution}°)...")
        
        count = 0
        for lat in range(-90, 91, resolution):
            for lon in range(-180, 181, resolution):
                observer_pos = (lat, lon, 0)  # Sea level
                dop = self.compute_dop(observer_pos, timestamp, mask_angle)
                
                lats.append(lat)
                lons.append(lon)
                dop_values.append(dop)
                
                count += 1
                if count % 100 == 0:
                    print(f"  Progress: {count}/{total_points} points ({100*count/total_points:.1f}%)")
        
        print(f"✓ Completed {count} points")
        return lats, lons, dop_values