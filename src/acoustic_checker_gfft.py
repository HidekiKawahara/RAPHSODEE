#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Acoustic Checker by Giant-FFT (PySide6 + PyQtGraph Desktop Implementation)
Ported from acouticCheckerByGFFTRefactored_v2_exportedR.m and constantQgaussBPFweight.m

MIT License

Copyright (c) 2026 Hideki Kawahara

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys
import time
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import Optional

import numpy as np
import scipy.fft
from scipy.special import erfinv
from scipy.interpolate import interp1d
import sounddevice as sd
import soundfile as sf

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QGridLayout,
    QVBoxLayout,
    QHBoxLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QComboBox,
    QRadioButton,
    QSlider,
    QDoubleSpinBox,
    QFileDialog,
    QMessageBox,
)
import pyqtgraph as pg


# ==============================================================================
# 1. DSP Core: Port of constantQgaussBPFweight.m (Log-Frequency Gaussian BPF)
# ==============================================================================
@dataclass
class ConstantQBPFWeight:
    nChannel: int
    filterWeight: np.ndarray      # (nChannel, nFreq) normalized weights
    filterWeightRaw: np.ndarray   # (nChannel, nFreq) raw weights
    quadfilterWeight: np.ndarray  # (nFreq, nChannel) sqrt(abs(raw)) for quadrature filtering
    fcSet: np.ndarray             # (nChannel,) center frequencies in Hz
    bandWidthInOctave: float
    stepInOctave: float
    elapsedTime: float


def constant_q_gauss_bpf_weight(
    frequency_axis: np.ndarray,
    bandwidth_in_octave: float,
    step_in_octave: float,
) -> ConstantQBPFWeight:
    """
    Exact NumPy port of constantQgaussBPFweight.m (log2-frequency Gaussian weights
    used by acouticCheckerByGFFTRefactored_v2_exportedR.m).
    """
    t0 = time.perf_counter()
    freq_ax = np.asarray(frequency_axis, dtype=np.float64).ravel()
    n_freq = freq_ax.size
    if n_freq < 2:
        raise ValueError("frequency_axis must contain at least two elements.")

    # Effective sampling upper limit consistent with MATLAB original:
    # fsEstimate = frequencyAxis(end) + frequencyAxis(2)
    fs_estimate = freq_ax[-1] + freq_ax[1]

    fc_low = 1000.0 * (2.0 ** -6)
    max_octave_index = np.log2((fs_estimate / 2.0) / fc_low)
    octave_indices = np.arange(0.0, max_octave_index + 1e-12, step_in_octave)
    fc_set = fc_low * (2.0 ** octave_indices)
    n_channel = fc_set.size

    # norminv(0.75) * 2 = sqrt(2) * erfinv(0.5) * 2
    base_sigma = np.sqrt(2.0) * erfinv(0.5) * 2.0
    sigma_log = bandwidth_in_octave / 2.0 / base_sigma

    # Build Gaussian weights on the log2-frequency axis and normalize per channel.
    # For f == 0 Hz, log2(0) is -inf so exp(-(-inf)^2) = 0.0.
    g_raw = np.zeros((n_freq, n_channel), dtype=np.float64)
    pos_mask = freq_ax > 0.0
    if np.any(pos_mask):
        log_freq_pos = np.log2(freq_ax[pos_mask])[:, None]       # (n_pos, 1)
        center_log = np.log2(fc_set)[None, :]                    # (1, n_channel)
        diff_log = log_freq_pos - center_log                     # (n_pos, n_channel)
        g_raw[pos_mask, :] = np.exp(-(diff_log ** 2) / (2.0 * (sigma_log ** 2)))

    col_sums = np.maximum(np.sum(g_raw, axis=0, keepdims=True), 1e-15)
    g_norm = g_raw / col_sums

    elapsed = time.perf_counter() - t0
    return ConstantQBPFWeight(
        nChannel=n_channel,
        filterWeight=g_norm.T.copy(),             # (nChannel, nFreq)
        filterWeightRaw=g_raw.T.copy(),           # (nChannel, nFreq)
        quadfilterWeight=np.sqrt(np.abs(g_raw)),  # (nFreq, nChannel)
        fcSet=fc_set,
        bandWidthInOctave=bandwidth_in_octave,
        stepInOctave=step_in_octave,
        elapsedTime=elapsed,
    )


# ==============================================================================
# 2. Main Desktop GUI & Real-Time Audio Engine
# ==============================================================================
class AcousticCheckerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Acoustic Checker by Giant-FFT (Python / PySide6)")
        self.resize(1220, 720)

        # Configure PyQtGraph global look
        pg.setConfigOptions(antialias=False, background="w", foreground="k")

        # --- State Properties (matching MATLAB App properties) ---
        self.sample_rate: int = 48000
        self.frame_length: int = 1024
        self.n_channel_in: int = 1
        self.n_channel_out: int = 2
        self.maximum_duration: int = 300  # seconds

        self.working_directory: Optional[Path] = None
        self.selected_device_id: Optional[int] = None
        self.selected_device_name: str = ""
        self.supported_rates: list[int] = []

        # Cyclic Audio Data
        self.cyclic_file_name: str = ""
        self.cyclic_file_path: Optional[Path] = None
        self.cyclic_data: Optional[np.ndarray] = None
        self.cyclic_fs: int = 48000

        # Streaming & Ring Buffer State
        self.stream: Optional[sd.Stream] = None
        self.streaming_enable: bool = False
        self.whole_source_buffer: Optional[np.ndarray] = None
        self.captured_buffer: Optional[np.ndarray] = None
        self.cycle_length: int = 0
        self.head_idx: int = 0
        self.play_idx: int = 0
        self.first_time: bool = True

        self.initial_source: Optional[np.ndarray] = None
        self.initial_captured: Optional[np.ndarray] = None

        # Analysis & Filter Bank State
        self.impulse_length: int = 48000
        self.data_length: int = 48000
        self.half_data_length: int = 24000
        self.n_fft: int = 48000
        self.fx_impulse: np.ndarray = np.arange(48000, dtype=np.float64)
        self.fxx_data: np.ndarray = np.arange(48000, dtype=np.float64)
        self.gf_str6: Optional[ConstantQBPFWeight] = None
        self.gf_str12: Optional[ConstantQBPFWeight] = None
        self.bpf_fc: float = 1000.0
        self.fc_1khz_filter_f: Optional[np.ndarray] = None
        self.fc_slider_filter_f: Optional[np.ndarray] = None

        # Cached latest spectra for instant radio-button / slider refresh
        self.last_lti_spec: Optional[np.ndarray] = None
        self.last_rtv_spec: Optional[np.ndarray] = None

        # Calibration State
        self.cal_initialized: bool = False
        self.ch1_cal_db: float = 0.0
        self.ch1_measured_spl: float = 94.0
        self.cal_signal: Optional[np.ndarray] = None
        self.cal_recorded_chunks: list[np.ndarray] = []
        self.cal_start_time: float = 0.0

        # Voice Memo State
        self.memo_start_time: float = 0.0
        self.memo_backup_msg: str = ""

        # Timers
        self.stream_timer = QTimer(self)
        self.stream_timer.setInterval(60)  # 60 ms polling period (same as MATLAB)
        self.stream_timer.timeout.connect(self.analysis_and_update_graphics_handler)

        self.memo_timer = QTimer(self)
        self.memo_timer.setInterval(100)
        self.memo_timer.timeout.connect(self._on_memo_timer_tick)

        self.cal_timer = QTimer(self)
        self.cal_timer.setInterval(100)
        self.cal_timer.timeout.connect(self._on_cal_timer_tick)

        # Build UI and initialize plots
        self._build_ui()
        self._initialize_graphics()
        self._initialize_working_directory()

    # ==========================================================================
    # UI Construction (3-Column Layout: Left, Center, Right)
    # ==========================================================================
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_grid = QGridLayout(central)
        main_grid.setContentsMargins(8, 8, 8, 8)
        main_grid.setSpacing(10)

        # ---------------- LEFT PANEL ----------------
        left_panel = QGroupBox("System & Stream Control")
        left_panel.setFixedWidth(285)
        left_layout = QVBoxLayout(left_panel)

        lbl_sys_title = QLabel("System Status")
        lbl_sys_title.setFont(QFont("Arial", 11, QFont.Bold))
        left_layout.addWidget(lbl_sys_title)

        self.status_label = QLabel("Click 'Set working folder' to begin.")
        self.status_label.setWordWrap(True)
        self.status_label.setMinimumHeight(55)
        self.status_label.setStyleSheet(
            "QLabel { background-color: #f4f6f8; border: 1px solid #ccc; padding: 6px; border-radius: 4px; }"
        )
        left_layout.addWidget(self.status_label)

        # IO Device Selection
        self.io_device_label = QLabel("IO device selection")
        self.io_device_label.setEnabled(False)
        left_layout.addWidget(self.io_device_label)

        self.io_device_dropdown = QComboBox()
        self.io_device_dropdown.setEnabled(False)
        self.io_device_dropdown.currentIndexChanged.connect(self.on_io_device_changed)
        left_layout.addWidget(self.io_device_dropdown)

        # Sampling Frequency Selection
        self.fs_label = QLabel("Sampling frequency")
        self.fs_label.setEnabled(False)
        left_layout.addWidget(self.fs_label)

        self.fs_dropdown = QComboBox()
        self.fs_dropdown.setEnabled(False)
        self.fs_dropdown.currentIndexChanged.connect(self.on_sampling_frequency_changed)
        left_layout.addWidget(self.fs_dropdown)

        # Input Channels
        self.in_ch_label = QLabel("Input Channels")
        self.in_ch_label.setEnabled(False)
        left_layout.addWidget(self.in_ch_label)

        self.in_ch_dropdown = QComboBox()
        self.in_ch_dropdown.addItems(["0", "1", "2"])
        self.in_ch_dropdown.setCurrentText("0")
        self.in_ch_dropdown.setEnabled(False)
        self.in_ch_dropdown.currentIndexChanged.connect(self.on_input_channel_changed)
        left_layout.addWidget(self.in_ch_dropdown)

        # Load Cyclic Data Button
        self.btn_load_cyclic = QPushButton("Load cyclic data")
        self.btn_load_cyclic.setEnabled(False)
        self.btn_load_cyclic.clicked.connect(self.on_load_cyclic_data_pushed)
        left_layout.addWidget(self.btn_load_cyclic)

        # Output Channel Group
        self.out_ch_group = QGroupBox("Output channel")
        self.out_ch_group.setEnabled(False)
        out_ch_layout = QHBoxLayout(self.out_ch_group)
        self.rb_out_l = QRadioButton("L")
        self.rb_out_r = QRadioButton("R")
        self.rb_out_lr = QRadioButton("L and R")
        self.rb_out_l.setChecked(True)
        out_ch_layout.addWidget(self.rb_out_l)
        out_ch_layout.addWidget(self.rb_out_r)
        out_ch_layout.addWidget(self.rb_out_lr)
        left_layout.addWidget(self.out_ch_group)

        # Spectrum Display Group
        self.spec_disp_group = QGroupBox("Spectrum display")
        self.spec_disp_group.setEnabled(False)
        spec_disp_layout = QVBoxLayout(self.spec_disp_group)
        self.rb_raw_dft = QRadioButton("Raw DFT")
        self.rb_oct_6th = QRadioButton("1/3 octave smoothed")
        self.rb_oct_12th = QRadioButton("1/12 octave smoothed")
        self.rb_raw_dft.setChecked(True)
        self.rb_raw_dft.toggled.connect(self._update_transfer_function_display)
        self.rb_oct_6th.toggled.connect(self._update_transfer_function_display)
        self.rb_oct_12th.toggled.connect(self._update_transfer_function_display)
        spec_disp_layout.addWidget(self.rb_raw_dft)
        spec_disp_layout.addWidget(self.rb_oct_6th)
        spec_disp_layout.addWidget(self.rb_oct_12th)
        left_layout.addWidget(self.spec_disp_group)

        # START / STOP Buttons
        self.btn_start = QPushButton("START")
        self.btn_start.setEnabled(False)
        self.btn_start.setMinimumHeight(32)
        self.btn_start.clicked.connect(self.on_start_pushed)
        left_layout.addWidget(self.btn_start)

        self.btn_stop = QPushButton("STOP")
        self.btn_stop.setEnabled(False)
        self.btn_stop.setMinimumHeight(32)
        self.btn_stop.clicked.connect(self.on_stop_pushed)
        left_layout.addWidget(self.btn_stop)

        left_layout.addStretch()

        self.btn_snapshot_gui = QPushButton("snapshot GUI")
        self.btn_snapshot_gui.setEnabled(False)
        self.btn_snapshot_gui.clicked.connect(self.on_snapshot_gui_pushed)
        left_layout.addWidget(self.btn_snapshot_gui)

        # ---------------- CENTER PANEL ----------------
        center_panel = QGroupBox("Workspace & Calibration")
        center_panel.setFixedWidth(185)
        center_layout = QVBoxLayout(center_panel)

        self.btn_set_folder = QPushButton("Set working folder")
        self.btn_set_folder.clicked.connect(self.on_set_working_folder_pushed)
        center_layout.addWidget(self.btn_set_folder)

        # Recording Lamp
        lamp_row = QHBoxLayout()
        self.lbl_rec_lamp = QLabel("Recording")
        self.lbl_rec_lamp.setEnabled(False)
        self.rec_lamp = QLabel()
        self.rec_lamp.setFixedSize(18, 18)
        self._set_lamp_active(False)
        lamp_row.addWidget(self.lbl_rec_lamp)
        lamp_row.addWidget(self.rec_lamp)
        center_layout.addLayout(lamp_row)

        self.btn_voice_memo = QPushButton("Voice memo 5s")
        self.btn_voice_memo.setEnabled(False)
        self.btn_voice_memo.clicked.connect(self.on_voice_memo_pushed)
        center_layout.addWidget(self.btn_voice_memo)

        # Calibration Ch-1 Panel
        self.cal_ch1_group = QGroupBox("Calibration Ch-1")
        self.cal_ch1_group.setEnabled(False)
        cal1_layout = QVBoxLayout(self.cal_ch1_group)

        self.btn_calibrate_1 = QPushButton("Calibrate")
        self.btn_calibrate_1.setEnabled(False)
        self.btn_calibrate_1.clicked.connect(self.on_calibrate_pushed)
        cal1_layout.addWidget(self.btn_calibrate_1)

        self.lbl_measured_spl = QLabel("Measured SPL")
        self.lbl_measured_spl.setEnabled(False)
        cal1_layout.addWidget(self.lbl_measured_spl)

        self.measured_spl_dropdown = QComboBox()
        self.measured_spl_dropdown.setEnabled(False)
        self.measured_spl_dropdown.currentIndexChanged.connect(self.on_measured_spl_changed)
        cal1_layout.addWidget(self.measured_spl_dropdown)

        self.lbl_cal_val = QLabel("CAL. Value: dB")
        self.lbl_cal_val.setEnabled(False)
        cal1_layout.addWidget(self.lbl_cal_val)

        self.cal_val_spin = QDoubleSpinBox()
        self.cal_val_spin.setRange(-200.0, 200.0)
        self.cal_val_spin.setDecimals(1)
        self.cal_val_spin.setEnabled(False)
        cal1_layout.addWidget(self.cal_val_spin)

        center_layout.addWidget(self.cal_ch1_group)
        center_layout.addStretch()

        # ---------------- RIGHT PANEL (Plots) ----------------
        right_panel = QGroupBox("Real-Time Acoustic Analysis")
        right_layout = QVBoxLayout(right_panel)

        # Top Plot: SpectrumUIAxes (Transfer Function Display)
        self.spectrum_plot = pg.PlotWidget(title="  ")
        self.spectrum_plot.setLogMode(x=True, y=False)
        self.spectrum_plot.showGrid(x=True, y=True, alpha=0.3)
        self.spectrum_plot.setLabel("left", "Level (rel. dB)")
        self.spectrum_plot.setLabel("bottom", "frequency (Hz)")
        self.spectrum_plot.setXRange(np.log10(20.0), np.log10(self.sample_rate / 2.0), padding=0.0)
        self.spectrum_plot.setYRange(-80.0, 20.0, padding=0.0)
        self._configure_transfer_function_x_ticks()
        right_layout.addWidget(self.spectrum_plot, stretch=5)

        # Slider Row: BPF fc
        slider_row = QHBoxLayout()
        self.lbl_bpf_slider = QLabel("BPF fc: 1000 Hz")
        self.lbl_bpf_slider.setMinimumWidth(115)
        self.bpf_slider = QSlider(Qt.Horizontal)
        self.bpf_slider.setMinimum(1)
        self.bpf_slider.setMaximum(100)
        self.bpf_slider.setValue(50)
        self.bpf_slider.setEnabled(False)
        self.bpf_slider.valueChanged.connect(self.on_bpf_slider_changed)
        slider_row.addWidget(self.lbl_bpf_slider)
        slider_row.addWidget(self.bpf_slider)
        right_layout.addLayout(slider_row)

        # Bottom Plots Row: BandSPLAxes (left) and bpfPowerAxes (right)
        bottom_plots_row = QHBoxLayout()

        self.band_spl_plot = pg.PlotWidget(title="  ")
        self.band_spl_plot.setLogMode(x=True, y=False)
        self.band_spl_plot.showGrid(x=True, y=True, alpha=0.3)
        self.band_spl_plot.setLabel("left", "SPL CAL.NO (dB)")
        self.band_spl_plot.setLabel("bottom", "frequency (Hz)")
        self.band_spl_plot.setXRange(np.log10(20.0), np.log10(self.sample_rate / 2.0), padding=0.0)
        self.band_spl_plot.setYRange(-110.0, 0.0, padding=0.0)
        bottom_plots_row.addWidget(self.band_spl_plot, stretch=1)

        self.power_resp_plot = pg.PlotWidget(title="Power response")
        self.power_resp_plot.showGrid(x=True, y=True, alpha=0.3)
        self.power_resp_plot.setLabel("left", "level (dB)")
        self.power_resp_plot.setLabel("bottom", "time (s)")
        self.power_resp_plot.setXRange(0.0, 1.0, padding=0.0)
        self.power_resp_plot.setYRange(-70.0, 0.0, padding=0.0)
        bottom_plots_row.addWidget(self.power_resp_plot, stretch=1)

        right_layout.addLayout(bottom_plots_row, stretch=4)

        # Snapshot Level Button
        snap_row = QHBoxLayout()
        self.btn_snapshot_level = QPushButton("Snapshot Level")
        self.btn_snapshot_level.setEnabled(False)
        self.btn_snapshot_level.clicked.connect(self.on_snapshot_level_pushed)
        snap_row.addWidget(self.btn_snapshot_level)
        snap_row.addStretch()
        right_layout.addLayout(snap_row)

        # Assemble 3 Panels
        main_grid.addWidget(left_panel, 0, 0)
        main_grid.addWidget(center_panel, 0, 1)
        main_grid.addWidget(right_panel, 0, 2)
        main_grid.setColumnStretch(2, 1)

    def _configure_transfer_function_x_ticks(self):
        """
        Labels only 100 Hz, 1000 Hz, and 10000 Hz on the transfer function display
        while keeping standard logarithmic minor grid lines.
        """
        major_ticks = [
            (2.0, "100 Hz"),
            (3.0, "1000 Hz"),
            (4.0, "10000 Hz"),
        ]
        minor_ticks = []
        for decade_exp in (1, 2, 3, 4):
            for mult in range(2, 10):
                val = (10 ** decade_exp) * mult
                if 20 <= val <= 100000:
                    minor_ticks.append((float(np.log10(val)), ""))
        bottom_axis = self.spectrum_plot.getAxis("bottom")
        bottom_axis.setTicks([major_ticks, minor_ticks])

    def _set_lamp_active(self, active: bool):
        if active:
            self.rec_lamp.setStyleSheet(
                "background-color: #2ecc71; border-radius: 9px; border: 1px solid #27ae60;"
            )
            self.lbl_rec_lamp.setEnabled(True)
        else:
            self.rec_lamp.setStyleSheet(
                "background-color: #bdc3c7; border-radius: 9px; border: 1px solid #95a5a6;"
            )
            self.lbl_rec_lamp.setEnabled(False)

    # ==========================================================================
    # Graphics Initialization & Helpers
    # ==========================================================================
    @staticmethod
    def _set_log_x_curve(curve: pg.PlotDataItem, fx: np.ndarray, y: np.ndarray):
        """
        Safely updates a pyqtgraph PlotDataItem on a log-X axis by masking f <= 0.
        """
        valid = fx > 0
        if np.any(valid):
            curve.setData(fx[valid], y[valid])

    def _initialize_graphics(self):
        fx = np.arange(1, self.sample_rate, dtype=np.float64)
        zeros_vec = np.zeros_like(fx)

        # Spectrum curves
        self.lti_curve = self.spectrum_plot.plot(
            fx, zeros_vec, pen=pg.mkPen(color=(0, 114, 189), width=2), name="LTI"
        )
        self.rtv_curve = self.spectrum_plot.plot(
            fx, zeros_vec - 30.0, pen=pg.mkPen(color=(217, 83, 25), width=2), name="RTV"
        )
        self.sdti_curve = self.spectrum_plot.plot(
            fx, zeros_vec - 20.0, pen=pg.mkPen(color=(237, 177, 32), width=2), name="SDTI"
        )
        self.spec_fc_line = pg.InfiniteLine(
            pos=np.log10(1000.0),
            angle=90,
            pen=pg.mkPen("g", width=3),
            label="1000",
            labelOpts={"position": 0.9, "color": (0, 140, 0)},
        )
        self.spectrum_plot.addItem(self.spec_fc_line)

        # BandSPL curves
        self.sig_lvl_curve = self.band_spl_plot.plot(
            fx, zeros_vec, pen=pg.mkPen(color=(0, 114, 189), width=2)
        )
        self.dist_lvl_curve = self.band_spl_plot.plot(
            fx, zeros_vec - 100.0, pen=pg.mkPen(color=(217, 83, 25), width=2)
        )
        self.ind_lvl_curve = self.band_spl_plot.plot(
            fx, zeros_vec - 100.0, pen=pg.mkPen(color=(237, 177, 32), width=2)
        )
        self.snap_lvl_curve = self.band_spl_plot.plot(
            fx, zeros_vec - 100.0, pen=pg.mkPen("g", width=2)
        )
        self.band_fc_line = pg.InfiniteLine(
            pos=np.log10(1000.0),
            angle=90,
            pen=pg.mkPen("g", width=3),
            label="1000",
            labelOpts={"position": 0.9, "color": (0, 140, 0)},
        )
        self.band_spl_plot.addItem(self.band_fc_line)

        # Power response curves
        tresp = np.arange(1, self.sample_rate + 1, dtype=np.float64) / self.sample_rate
        self.presp_1k_curve = self.power_resp_plot.plot(
            tresp, np.zeros_like(tresp), pen=pg.mkPen(color=(0, 114, 189), width=2)
        )
        self.presp_fc_curve = self.power_resp_plot.plot(
            tresp, np.zeros_like(tresp), pen=pg.mkPen("g", width=2)
        )

    def _update_graphics_axes(self):
        fs = self.sample_rate
        fftl = self.impulse_length
        fx = np.arange(fftl, dtype=np.float64) / fftl * fs
        self.fx_impulse = fx
        self.fxx_data = np.arange(self.data_length, dtype=np.float64) / self.data_length * fs

        self.spectrum_plot.setXRange(np.log10(20.0), np.log10(fs / 2.0), padding=0.0)
        self.band_spl_plot.setXRange(np.log10(20.0), np.log10(fs / 2.0), padding=0.0)
        self._configure_transfer_function_x_ticks()

        self._set_log_x_curve(self.lti_curve, fx, np.zeros(fftl))
        self._set_log_x_curve(self.rtv_curve, fx, np.zeros(fftl) - 30.0)
        self._set_log_x_curve(self.sdti_curve, fx, np.zeros(fftl) - 120.0)

        self._set_log_x_curve(self.sig_lvl_curve, fx, np.zeros(fftl))
        self._set_log_x_curve(self.dist_lvl_curve, fx, np.zeros(fftl) - 200.0)
        self._set_log_x_curve(self.ind_lvl_curve, fx, np.zeros(fftl) - 200.0)
        self._set_log_x_curve(self.snap_lvl_curve, fx, np.zeros(fftl) - 200.0)

        # Configure BPF slider bounds with +/-11 guard margin (matching MATLAB 1-based [12, nBands-11])
        if self.gf_str12 is not None:
            n_bands = self.gf_str12.nChannel
            fc_set = self.gf_str12.fcSet
            idx_1khz = int(np.argmin(np.abs(fc_set - 1000.0)))
            min_idx = min(11, n_bands - 1)
            max_idx = max(min_idx, n_bands - 12)
            self.bpf_slider.blockSignals(True)
            self.bpf_slider.setMinimum(min_idx)
            self.bpf_slider.setMaximum(max_idx)
            clamped_val = min(max(idx_1khz, min_idx), max_idx)
            self.bpf_slider.setValue(clamped_val)
            self.bpf_slider.blockSignals(False)
            self.bpf_slider.setEnabled(True)

    # ==========================================================================
    # Workflow & Setup Callbacks
    # ==========================================================================
    def _initialize_working_directory(self):
        self.status_label.setText("Click 'Set working folder' to begin.")
        self.btn_set_folder.setEnabled(True)
        self.io_device_dropdown.setEnabled(False)
        self.io_device_label.setEnabled(False)

    def on_set_working_folder_pushed(self):
        sel_dir = QFileDialog.getExistingDirectory(self, "Select Working Directory")
        if not sel_dir:
            return
        tstr = datetime.now().strftime("%Y%m%d%H%M%S")
        work_dir = Path(sel_dir) / f"acTest{tstr}"
        work_dir.mkdir(parents=True, exist_ok=True)
        self.working_directory = work_dir

        self.btn_set_folder.setEnabled(False)
        self.btn_snapshot_gui.setEnabled(True)
        self.populate_audio_devices()

    def populate_audio_devices(self):
        self.fs_dropdown.setEnabled(False)
        self.in_ch_dropdown.setEnabled(False)
        self.btn_load_cyclic.setEnabled(False)
        self.btn_start.setEnabled(False)
        self.btn_stop.setEnabled(False)
        self.streaming_enable = False

        try:
            devices = sd.query_devices()
        except Exception as e:
            QMessageBox.critical(self, "Audio Error", f"Cannot query audio devices: {e}")
            return

        # Filter full-duplex devices (both input and output channels > 0)
        duplex_items = []
        for idx, dev in enumerate(devices):
            if dev.get("max_input_channels", 0) > 0 and dev.get("max_output_channels", 0) > 0:
                duplex_items.append((idx, f"{idx}: {dev['name']}"))

        self.io_device_dropdown.blockSignals(True)
        self.io_device_dropdown.clear()
        if not duplex_items:
            self.io_device_dropdown.addItem("<No full-duplex devices found>", -1)
            self.io_device_dropdown.setEnabled(False)
            self.status_label.setText("No full-duplex audio devices available")
            self.io_device_dropdown.blockSignals(False)
            return

        for idx, label in duplex_items:
            self.io_device_dropdown.addItem(label, idx)
        self.io_device_dropdown.addItem(" ---- ", -1)
        self.io_device_dropdown.setCurrentIndex(self.io_device_dropdown.count() - 1)
        self.io_device_dropdown.blockSignals(False)

        self.io_device_dropdown.setEnabled(True)
        self.io_device_label.setEnabled(True)
        self.status_label.setText("Select IO device")

    def on_io_device_changed(self, index: int):
        dev_id = self.io_device_dropdown.itemData(index)
        if dev_id is None or dev_id < 0:
            return

        self.selected_device_id = int(dev_id)
        dev_info = sd.query_devices(self.selected_device_id)
        self.selected_device_name = dev_info["name"]
        self.status_label.setText(f"Selected: {self.selected_device_name}")

        candidate_rates = [8000, 11025, 16000, 22050, 32000, 44100, 48000, 88200, 96000, 192000]
        supported = []
        for r in candidate_rates:
            try:
                sd.check_input_settings(device=self.selected_device_id, samplerate=r, channels=1)
                sd.check_output_settings(device=self.selected_device_id, samplerate=r, channels=1)
                supported.append(r)
            except Exception:
                pass

        if not supported:
            QMessageBox.warning(
                self, "Device Error", f"No supported full-duplex rates found for {self.selected_device_name}"
            )
            return

        self.supported_rates = supported
        self.fs_dropdown.blockSignals(True)
        self.fs_dropdown.clear()
        for r in supported:
            self.fs_dropdown.addItem(str(r), r)
        self.fs_dropdown.addItem(" ---- ", -1)
        self.fs_dropdown.setCurrentIndex(self.fs_dropdown.count() - 1)
        self.fs_dropdown.blockSignals(False)

        self.fs_dropdown.setEnabled(True)
        self.fs_label.setEnabled(True)
        self.io_device_dropdown.setEnabled(False)

        max_out_ch = dev_info.get("max_output_channels", 1)
        if max_out_ch > 1:
            self.n_channel_out = 2
            self.out_ch_group.setEnabled(True)
        else:
            self.n_channel_out = 1
            self.out_ch_group.setEnabled(False)

        self.status_label.setText("Select Sampling frequency")

    def on_sampling_frequency_changed(self, index: int):
        r = self.fs_dropdown.itemData(index)
        if r is None or r < 0:
            return

        self.sample_rate = int(r)
        self.status_label.setText(f"Sampling frequency is {self.sample_rate} Hz")
        self.fs_dropdown.setEnabled(False)

        dev_info = sd.query_devices(self.selected_device_id)
        max_rec_ch = dev_info.get("max_input_channels", 1)
        if max_rec_ch > 1:
            self.in_ch_dropdown.setEnabled(True)
            self.in_ch_label.setEnabled(True)
            self.status_label.setText("Select Number of Input Channels")
        else:
            self.in_ch_dropdown.setCurrentText("1")
            self.in_ch_dropdown.setEnabled(False)
            self.on_input_channel_changed(self.in_ch_dropdown.currentIndex())

    def on_input_channel_changed(self, index: int):
        val_str = self.in_ch_dropdown.currentText()
        n_ch = int(val_str) if val_str.isdigit() else 0
        if n_ch > 0:
            self.n_channel_in = n_ch
            self.in_ch_dropdown.setEnabled(False)
            self.status_label.setText("Ready to load cyclic data")
            self.btn_load_cyclic.setEnabled(True)
            self.spec_disp_group.setEnabled(True)
            self.btn_voice_memo.setEnabled(True)

    # ==========================================================================
    # Load Cyclic Data & Design Filter Banks
    # ==========================================================================
    def on_load_cyclic_data_pushed(self):
        fpath, _ = QFileDialog.getOpenFileName(
            self, "Select Cyclic Audio File", "", "Audio Files (*.wav *.aif *.aiff *.WAV *.AIF)"
        )
        if not fpath:
            return

        try:
            x, fs = sf.read(fpath, dtype="float64", always_2d=True)
        except Exception as e:
            QMessageBox.critical(self, "File Error", f"Failed to read audio file: {e}")
            return

        if x.size == 0:
            QMessageBox.warning(self, "File Error", "Selected audio file is empty.")
            return

        # Enforce mono column vector
        x_mono = x[:, 0].copy()
        self.cyclic_file_path = Path(fpath)
        self.cyclic_file_name = self.cyclic_file_path.name
        self.cyclic_data = x_mono
        self.cyclic_fs = int(fs)

        if self.sample_rate != self.cyclic_fs:
            if self.cyclic_fs not in self.supported_rates:
                self.status_label.setText("This sampling rate is not supported. Select another file.")
                return
            self.sample_rate = self.cyclic_fs
            idx = self.fs_dropdown.findData(self.sample_rate)
            if idx >= 0:
                self.fs_dropdown.blockSignals(True)
                self.fs_dropdown.setCurrentIndex(idx)
                self.fs_dropdown.blockSignals(False)

        self.spectrum_plot.setTitle(f"File: {self.cyclic_file_name}")
        self.impulse_length = min(self.sample_rate, int(np.floor(x_mono.size / 3)))
        self.data_length = x_mono.size
        self.half_data_length = int(np.floor(x_mono.size / 2))
        self.fx_impulse = (
            np.arange(self.impulse_length, dtype=np.float64) / self.impulse_length * self.sample_rate
        )

        # Build Constant-Q Log-Gaussian BPF weights (1/3 octave and 1/12 octave)
        # matching MATLAB's constantQgaussBPFweight(frequencyAxisRow, 1/3, 1/24)
        self.gf_str6 = constant_q_gauss_bpf_weight(self.fx_impulse, 1.0 / 3.0, 1.0 / 24.0)
        self.gf_str12 = constant_q_gauss_bpf_weight(self.fx_impulse, 1.0 / 12.0, 1.0 / 24.0)

        self.bpf_fc = 1000.0
        tmp_bpf = self._design_fc_octave_set(self.bpf_fc)
        self.fc_1khz_filter_f = tmp_bpf
        self.fc_slider_filter_f = tmp_bpf

        self.last_lti_spec = None
        self.last_rtv_spec = None

        self._update_graphics_axes()
        self.status_label.setText("Ready to START")
        self.btn_start.setEnabled(True)

    def _design_fc_octave_set(self, fc: float) -> np.ndarray:
        fc_set = self.gf_str12.fcSet
        n_bands = fc_set.size
        idx_fc = int(np.argmin(np.abs(fc - fc_set)))
        band_indices = np.clip(idx_fc + np.arange(-11, 12), 0, n_bands - 1)
        return self.gf_str12.quadfilterWeight[:, band_indices].copy()

    def on_bpf_slider_changed(self, value: int):
        if self.gf_str12 is None:
            return
        fc_set = self.gf_str12.fcSet
        idx = int(np.clip(value, 0, fc_set.size - 1))
        fc = float(fc_set[idx])
        self.bpf_fc = fc
        self.lbl_bpf_slider.setText(f"BPF fc: {fc:6.0f} Hz")

        log_fc = np.log10(max(fc, 1.0))
        self.spec_fc_line.setPos(log_fc)
        self.spec_fc_line.label.setFormat(f"{fc:6.0f}")
        self.band_fc_line.setPos(log_fc)
        self.band_fc_line.label.setFormat(f"{fc:6.0f}")

        self.fc_slider_filter_f = self._design_fc_octave_set(fc)

        # If we already have a captured LTI spectrum, update green power response immediately
        if self.last_lti_spec is not None:
            eps_spec = 1e-12
            fftl = self.impulse_length
            tx = np.arange(1, fftl + 1, dtype=np.float64) / self.sample_rate
            resp_fc = scipy.fft.ifft(self.fc_slider_filter_f * self.last_lti_spec[:, None], axis=0)
            p_resp_fc = 10.0 * np.log10(np.maximum(np.sum(np.abs(resp_fc) ** 2, axis=1), eps_spec))
            p_resp_fc -= np.max(p_resp_fc)
            self.presp_fc_curve.setData(tx[::10], p_resp_fc[::10])

    def _update_transfer_function_display(self):
        """
        Updates the top Spectrum plot from the latest LTI and RTV spectra
        according to the currently selected radio button (Raw DFT, 1/3 oct, 1/12 oct).
        """
        if self.last_lti_spec is None or self.last_rtv_spec is None:
            return
        if self.gf_str6 is None or self.gf_str12 is None:
            return

        eps_spec = 1e-12
        lti_spec = self.last_lti_spec
        rtv_spec = self.last_rtv_spec
        lti_pw = np.abs(lti_spec) ** 2
        rtv_pw = np.abs(rtv_spec) ** 2

        if self.rb_raw_dft.isChecked():
            freq_axis = self.fx_impulse
            lti_gain = 20.0 * np.log10(np.maximum(np.abs(lti_spec), eps_spec))
            rtv_gain = 20.0 * np.log10(np.maximum(np.abs(rtv_spec), eps_spec))
        elif self.rb_oct_6th.isChecked():
            freq_axis = self.gf_str6.fcSet
            lti_gain = 10.0 * np.log10(np.maximum(self.gf_str6.filterWeight @ lti_pw, eps_spec))
            rtv_gain = 10.0 * np.log10(np.maximum(self.gf_str6.filterWeight @ rtv_pw, eps_spec))
        else:
            freq_axis = self.gf_str12.fcSet
            lti_gain = 10.0 * np.log10(np.maximum(self.gf_str12.filterWeight @ lti_pw, eps_spec))
            rtv_gain = 10.0 * np.log10(np.maximum(self.gf_str12.filterWeight @ rtv_pw, eps_spec))

        self._set_log_x_curve(self.lti_curve, freq_axis, lti_gain)
        self._set_log_x_curve(self.rtv_curve, freq_axis, rtv_gain)

    # ==========================================================================
    # Full-Duplex Streaming & Giant-FFT Real-Time Analysis
    # ==========================================================================
    def _prepare_output_cycle_buffer(self, unit_mono: np.ndarray) -> np.ndarray:
        if self.n_channel_out == 2:
            buf = np.column_stack([unit_mono, unit_mono])
            if self.rb_out_l.isChecked():
                buf[:, 1] = 0.0
            elif self.rb_out_r.isChecked():
                buf[:, 0] = 0.0
            return buf
        return unit_mono[:, None]

    def on_start_pushed(self):
        if self.cyclic_data is None:
            return

        self.btn_voice_memo.setEnabled(False)
        fs = self.sample_rate
        tmp_buffer = self._prepare_output_cycle_buffer(self.cyclic_data)
        self.cycle_length = tmp_buffer.shape[0]
        self.n_fft = self.cycle_length

        n_cycles = int(np.floor(self.maximum_duration * fs / self.cycle_length))
        if n_cycles < 1:
            n_cycles = 1

        self.whole_source_buffer = np.tile(tmp_buffer, (n_cycles, 1))
        total_samples = self.whole_source_buffer.shape[0]
        self.captured_buffer = np.zeros((total_samples, self.n_channel_in), dtype=np.float64)

        self.head_idx = 0
        self.play_idx = 0
        self.first_time = True
        self.streaming_enable = True

        def audio_callback(indata, outdata, frames, time_info, status):
            rem = total_samples - self.play_idx
            if rem <= 0:
                outdata.fill(0.0)
                raise sd.CallbackStop()
            n_copy = min(frames, rem)
            outdata[:n_copy, :] = self.whole_source_buffer[self.play_idx : self.play_idx + n_copy, :]
            if n_copy < frames:
                outdata[n_copy:, :].fill(0.0)

            rec_rem = total_samples - self.head_idx
            n_rec = min(frames, rec_rem)
            if n_rec > 0:
                self.captured_buffer[self.head_idx : self.head_idx + n_rec, :] = indata[:n_rec, :]
                self.head_idx += n_rec
            self.play_idx += n_copy

        try:
            self.stream = sd.Stream(
                samplerate=fs,
                blocksize=self.frame_length,
                device=(self.selected_device_id, self.selected_device_id),
                channels=(self.n_channel_in, self.n_channel_out),
                dtype="float32",
                callback=audio_callback,
            )
            self.stream.start()
        except Exception as e:
            QMessageBox.critical(self, "Audio Error", f"Cannot start audio stream: {e}")
            self.btn_voice_memo.setEnabled(True)
            return

        self.status_label.setText("Streaming")
        self.btn_start.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.btn_load_cyclic.setEnabled(False)
        self.out_ch_group.setEnabled(False)
        self.btn_snapshot_level.setEnabled(True)
        self.stream_timer.start()

    def analysis_and_update_graphics_handler(self):
        if self.stream is None or not self.stream.active:
            self.stream_timer.stop()
            self._stop_housekeeping()
            return

        curr_head = self.head_idx
        cycle_len = self.cycle_length
        total_len = self.whole_source_buffer.shape[0]

        if curr_head + cycle_len > total_len:
            self.on_stop_pushed()
            return

        if curr_head < cycle_len:
            return

        fs = self.sample_rate
        fxx = self.fxx_data
        fx_impulse = self.fx_impulse
        fc_set = self.gf_str6.fcSet
        fc_set_l = fc_set * (2.0 ** (-0.5))
        fc_set_h = fc_set * (2.0 ** 0.5)

        source_id = 1 if (self.rb_out_r.isChecked() and self.n_channel_out == 2) else 0
        start_idx = curr_head - cycle_len
        source_seg = self.whole_source_buffer[start_idx:curr_head, source_id].copy()
        captured_raw = self.captured_buffer[start_idx:curr_head, 0].copy()

        if self.first_time:
            captured_seg = source_seg.copy()
            self.initial_captured = captured_raw.copy()
            self.initial_source = source_seg.copy()
            self.first_time = False
        else:
            captured_seg = captured_raw

        eps_spec = 1e-12
        n_fft = self.n_fft

        captured_spectrum = scipy.fft.fft(captured_seg, n=n_fft)
        captured_power = (np.abs(captured_spectrum) ** 2) / (n_fft ** 2)

        cum_cap_pw = np.cumsum(captured_power)
        interp_cap = interp1d(fxx, cum_cap_pw, kind="linear", fill_value="extrapolate")
        octave_band_level = np.maximum(interp_cap(fc_set_h) - interp_cap(fc_set_l), eps_spec)

        transfer_fn = captured_spectrum / scipy.fft.fft(source_seg + eps_spec, n=n_fft)
        impulse_resp = np.real(scipy.fft.ifft(transfer_fn, n=n_fft))

        fftl = self.impulse_length
        n_imp = impulse_resp.size
        idx_peak = int(np.argmax(np.abs(impulse_resp)))
        bias_l = self.half_data_length

        # Modular indexing matching MATLAB:
        # indIndices = mod((0:fftl-1) + idxPeak + round(fs*0.015) - 1, Nimp) + 1
        # rtvIndices = mod((0:fftl-1) + biasL - 1, Nimp) + 1
        ind_indices = (np.arange(fftl) + idx_peak + int(np.floor(fs * 0.015 + 0.5))) % n_imp
        rtv_indices = (np.arange(fftl) + bias_l) % n_imp

        lti_spec = scipy.fft.fft(impulse_resp[:fftl])
        indirect_spec = scipy.fft.fft(impulse_resp[ind_indices])
        rtv_spec = scipy.fft.fft(impulse_resp[rtv_indices])

        self.last_lti_spec = lti_spec
        self.last_rtv_spec = rtv_spec

        lti_pw = np.abs(lti_spec) ** 2
        indirect_pw = np.abs(indirect_spec) ** 2
        rtv_pw = np.abs(rtv_spec) ** 2

        interp_lti = interp1d(fx_impulse, np.cumsum(lti_pw), kind="linear", fill_value="extrapolate")
        interp_rtv = interp1d(fx_impulse, np.cumsum(rtv_pw), kind="linear", fill_value="extrapolate")
        interp_ind = interp1d(fx_impulse, np.cumsum(indirect_pw), kind="linear", fill_value="extrapolate")

        oct_lti_lvl = np.maximum(interp_lti(fc_set_h) - interp_lti(fc_set_l), eps_spec)
        oct_rtv_lvl = np.maximum(interp_rtv(fc_set_h) - interp_rtv(fc_set_l), eps_spec)
        oct_ind_lvl = np.maximum(interp_ind(fc_set_h) - interp_ind(fc_set_l), eps_spec)

        rel_rtv_db = 10.0 * np.log10(oct_rtv_lvl / oct_lti_lvl)
        rel_ind_db = 10.0 * np.log10(oct_ind_lvl / oct_lti_lvl)

        # Update top Spectrum plot
        self._update_transfer_function_display()

        # Update bottom-left Band SPL plot
        oct_band_db = 10.0 * np.log10(octave_band_level) + self.ch1_cal_db + 3.0
        self._set_log_x_curve(self.sig_lvl_curve, fc_set, oct_band_db)
        self._set_log_x_curve(self.dist_lvl_curve, fc_set, oct_band_db + rel_rtv_db)
        self._set_log_x_curve(self.ind_lvl_curve, fc_set, oct_band_db + rel_ind_db)

        # Update bottom-right Power decay plot
        tx = np.arange(1, fftl + 1, dtype=np.float64) / fs
        resp_1k = scipy.fft.ifft(self.fc_1khz_filter_f * lti_spec[:, None], axis=0)
        p_resp_1k = 10.0 * np.log10(np.maximum(np.sum(np.abs(resp_1k) ** 2, axis=1), eps_spec))
        p_resp_1k -= np.max(p_resp_1k)
        self.presp_1k_curve.setData(tx[::10], p_resp_1k[::10])

        resp_fc = scipy.fft.ifft(self.fc_slider_filter_f * lti_spec[:, None], axis=0)
        p_resp_fc = 10.0 * np.log10(np.maximum(np.sum(np.abs(resp_fc) ** 2, axis=1), eps_spec))
        p_resp_fc -= np.max(p_resp_fc)
        self.presp_fc_curve.setData(tx[::10], p_resp_fc[::10])

    def on_stop_pushed(self):
        self.streaming_enable = False
        self.stream_timer.stop()
        self._close_audio_stream()
        self.btn_snapshot_level.setEnabled(False)

        if not self.cal_initialized:
            self._enable_calibration()

        self._stop_housekeeping()

        if self.working_directory is None:
            return

        tstr = datetime.now().strftime("%Y%m%d%H%M%S")

        # Save initial cycle
        if self.initial_source is not None and self.initial_captured is not None:
            sf.write(
                str(self.working_directory / f"initialSourceAT{tstr}r.wav"),
                self.initial_source,
                self.sample_rate,
                subtype="PCM_24",
            )
            sf.write(
                str(self.working_directory / f"initialCapturedAT{tstr}r.wav"),
                self.initial_captured,
                self.sample_rate,
                subtype="PCM_24",
            )

        # Save full accumulated stream
        if self.head_idx > 0 and self.whole_source_buffer is not None and self.captured_buffer is not None:
            end_idx = min(self.head_idx, self.whole_source_buffer.shape[0])
            sf.write(
                str(self.working_directory / f"sourceAT{tstr}r.wav"),
                self.whole_source_buffer[:end_idx, :],
                self.sample_rate,
                subtype="PCM_24",
            )
            sf.write(
                str(self.working_directory / f"capturedAT{tstr}r.wav"),
                self.captured_buffer[:end_idx, :],
                self.sample_rate,
                subtype="PCM_24",
            )

    def _stop_housekeeping(self):
        self._close_audio_stream()
        self.status_label.setText("Stopped: START or Select another file.")
        self.btn_stop.setEnabled(False)
        self.btn_start.setEnabled(True)
        self.btn_load_cyclic.setEnabled(True)
        self.btn_voice_memo.setEnabled(True)
        if self.n_channel_out > 1:
            self.out_ch_group.setEnabled(True)

    def _close_audio_stream(self):
        if self.stream is not None:
            try:
                self.stream.stop()
                self.stream.close()
            except Exception:
                pass
            self.stream = None

    # ==========================================================================
    # Calibration, Voice Memo, and Snapshots
    # ==========================================================================
    def _enable_calibration(self):
        self.cal_ch1_group.setEnabled(True)
        self.btn_calibrate_1.setEnabled(True)
        self.cal_val_spin.setEnabled(False)
        self.lbl_cal_val.setEnabled(False)
        self.cal_initialized = True

    def _generate_calibration_signal(self) -> np.ndarray:
        fs = self.sample_rate
        tt = np.arange(1, fs + 1, dtype=np.float64)[:, None] / fs
        fc = 1000.0
        fc_list = np.arange(round(fc * (2.0 ** -0.5)), round(fc * (2.0 ** 0.5)) + 1, dtype=np.float64)
        phases = np.random.rand(1, fc_list.size)
        x = np.sum(np.sin(2.0 * np.pi * (tt @ fc_list[None, :] + phases)), axis=1)
        unit_source = (x / np.std(x, ddof=1)) * 0.1
        return self._prepare_output_cycle_buffer(unit_source)

    def on_calibrate_pushed(self):
        spl_levels = [0] + list(range(95, 64, -1))
        self.measured_spl_dropdown.blockSignals(True)
        self.measured_spl_dropdown.clear()
        for lvl in spl_levels:
            self.measured_spl_dropdown.addItem(str(lvl), float(lvl))
        self.measured_spl_dropdown.setCurrentIndex(0)
        self.measured_spl_dropdown.blockSignals(False)

        self.lbl_measured_spl.setEnabled(True)
        unit_cal = self._generate_calibration_signal()
        cycle_len = unit_cal.shape[0]
        self.cal_recorded_chunks = []
        play_ptr = [0]

        def cal_callback(indata, outdata, frames, time_info, status):
            for i in range(frames):
                outdata[i, :] = unit_cal[play_ptr[0], :]
                play_ptr[0] = (play_ptr[0] + 1) % cycle_len
            self.cal_recorded_chunks.append(indata[:, 0].copy())

        try:
            self.stream = sd.Stream(
                samplerate=self.sample_rate,
                blocksize=self.frame_length,
                device=(self.selected_device_id, self.selected_device_id),
                channels=(self.n_channel_in, self.n_channel_out),
                dtype="float32",
                callback=cal_callback,
            )
            self.stream.start()
        except Exception as e:
            QMessageBox.critical(self, "Audio Error", f"Cannot start calibration stream: {e}")
            return

        self.memo_backup_msg = self.status_label.text()
        self.cal_start_time = time.perf_counter()
        self._set_lamp_active(True)
        self.cal_timer.start()

    def _on_cal_timer_tick(self):
        elapsed = time.perf_counter() - self.cal_start_time
        if elapsed > 3.1 and not self.measured_spl_dropdown.isEnabled():
            self.measured_spl_dropdown.setEnabled(True)
            self.status_label.setText("Select Measured SPL from dropdown to finish calibration.")

    def on_measured_spl_changed(self, index: int):
        val = self.measured_spl_dropdown.itemData(index)
        if val is None:
            return
        self.ch1_measured_spl = float(val)
        self.cal_timer.stop()
        self._close_audio_stream()
        self._set_lamp_active(False)

        if self.cal_recorded_chunks:
            y = np.concatenate(self.cal_recorded_chunks)
            tail_len = min(y.size, 3 * self.sample_rate)
            if tail_len > 0:
                observed_std_db = 20.0 * np.log10(max(float(np.std(y[-tail_len:], ddof=1)), 1e-12))
                cal_db = self.ch1_measured_spl - observed_std_db
                self.ch1_cal_db = cal_db
                self.cal_val_spin.setValue(cal_db)

        self.band_spl_plot.setYRange(0.0, 100.0, padding=0.0)
        self.band_spl_plot.setLabel("left", "CAL. SPL (dB)")
        self.btn_calibrate_1.setEnabled(False)
        self.measured_spl_dropdown.setEnabled(False)
        self.lbl_measured_spl.setEnabled(False)
        self.status_label.setText(self.memo_backup_msg)

    def on_voice_memo_pushed(self):
        fs = self.sample_rate
        n_frames = int(fs * 5.3)
        self.voice_memo_buf = np.zeros((n_frames, self.n_channel_in), dtype=np.float32)
        rec_ptr = [0]

        def memo_callback(indata, outdata, frames, time_info, status):
            outdata.fill(0.0)
            rem = n_frames - rec_ptr[0]
            if rem <= 0:
                raise sd.CallbackStop()
            n_copy = min(frames, rem)
            self.voice_memo_buf[rec_ptr[0] : rec_ptr[0] + n_copy, :] = indata[:n_copy, :]
            rec_ptr[0] += n_copy

        try:
            self.stream = sd.Stream(
                samplerate=fs,
                blocksize=self.frame_length,
                device=(self.selected_device_id, self.selected_device_id),
                channels=(self.n_channel_in, self.n_channel_out),
                dtype="float32",
                callback=memo_callback,
            )
            self.stream.start()
        except Exception as e:
            QMessageBox.critical(self, "Audio Error", f"Cannot start voice memo: {e}")
            return

        self.memo_backup_msg = self.status_label.text()
        self.memo_start_time = time.perf_counter()
        self._set_lamp_active(True)
        self.btn_voice_memo.setEnabled(False)
        self.memo_timer.start()

    def _on_memo_timer_tick(self):
        elapsed = time.perf_counter() - self.memo_start_time
        if elapsed <= 5.3:
            time_left = max(0.0, round((5.0 - elapsed) * 10.0) / 10.0)
            self.status_label.setText(f"Voice memo: {time_left:.1f} (s) left")
        else:
            self.memo_timer.stop()
            self._close_audio_stream()
            self._set_lamp_active(False)
            self.status_label.setText(self.memo_backup_msg)
            self.btn_voice_memo.setEnabled(True)

            if self.working_directory is not None:
                y = self.voice_memo_buf[: self.sample_rate * 5, :]
                max_amp = float(np.max(np.abs(y)))
                y_norm = (y / max_amp) * 0.5 if max_amp > 0 else y
                tstr = datetime.now().strftime("%Y%m%d%H%M%S")
                out_path = self.working_directory / f"voiceMemo{tstr}r.wav"
                sf.write(str(out_path), y_norm, self.sample_rate, subtype="PCM_24")

    def on_snapshot_level_pushed(self):
        x_data, y_data = self.sig_lvl_curve.getData()
        if x_data is not None and y_data is not None:
            self.snap_lvl_curve.setData(10.0 ** x_data, y_data)

    def on_snapshot_gui_pushed(self):
        if self.working_directory is None:
            return
        pixmap = self.grab()
        tstr = datetime.now().strftime("%Y%m%d%H%M%S")
        out_path = self.working_directory / f"GUIsnapAT{tstr}.png"
        pixmap.save(str(out_path), "PNG")

    def closeEvent(self, event):
        self.stream_timer.stop()
        self.memo_timer.stop()
        self.cal_timer.stop()
        self._close_audio_stream()
        event.accept()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = AcousticCheckerApp()
    window.show()
    sys.exit(app.exec())