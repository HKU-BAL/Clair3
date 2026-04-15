"""
TensorConfig - Extensible tensor configuration for Clair3 training.

Provides a centralized, config-driven way to manage tensor dimensions,
channels, and extra matrices. Supports:
  - Adding new channels (e.g., dwell_time, signal features)
  - Changing matrix depth (e.g., 89 → 128 for higher coverage)
  - Registering extra matrices alongside the main tensor (e.g., phasing_matrix)

Default values are read from shared/param_f.py and shared/param_p.py
to maintain backward compatibility.
"""

import numpy as np


class TensorComponent:
    """A single tensor component in the training data schema."""

    def __init__(self, name, shape, dtype='int8', optional=False):
        """
        Args:
            name: Dataset name in HDF5 (e.g., 'position_matrix', 'phasing_matrix')
            shape: Shape per sample, excluding batch dimension (e.g., [89, 33, 8])
            dtype: numpy dtype string (e.g., 'int8', 'int32', 'float16')
            optional: If True, this component may be absent in some files
        """
        self.name = name
        self.shape = list(shape) if shape else []
        self.dtype = dtype
        self.optional = optional

    def __repr__(self):
        return "TensorComponent({!r}, shape={}, dtype={!r}, optional={})".format(
            self.name, self.shape, self.dtype, self.optional)


class TensorConfig:
    """
    Extensible tensor configuration.

    Usage:
        # Standard full-alignment for ONT
        config = TensorConfig('full_alignment', 'ont')
        # → tensor_shape = [89, 33, 8]

        # With dwell time
        config = TensorConfig('full_alignment', 'ont', extra_channels=['dwell_time'])
        # → tensor_shape = [89, 33, 9]

        # Custom depth
        config = TensorConfig('full_alignment', 'ont', custom_depth=128)
        # → tensor_shape = [128, 33, 8]

        # Add extra matrix
        config.add_component('phasing_matrix', shape=[89, 64], dtype='int8', optional=True)
    """

    # Default values from param_f.py and param_p.py
    MATRIX_DEPTH_DICT = {'ont': 89, 'hifi': 55, 'ilmn': 55}
    FLANKING_BASE_NUM = 16
    NO_OF_POSITIONS = 33

    # Pileup: 18 channels (9 forward + 9 reverse)
    PILEUP_CHANNELS = ('A', 'C', 'G', 'T', 'I', 'I1', 'D', 'D1', '*',
                       'a', 'c', 'g', 't', 'i', 'i1', 'd', 'd1', '#')

    # Full-alignment: 8 base channels
    FA_CHANNELS = ('reference_base', 'alternative_base', 'mapping_quality',
                   'base_quality', 'strand_info', 'variant_type',
                   'insert_base', 'phasing_info')

    # Label: 21 GT + 3 zygosity + 33 len1 + 33 len2 = 90
    LABEL_SHAPE = [21, 3, 33, 33]
    LABEL_SIZE = sum(LABEL_SHAPE)

    def __init__(self, mode='full_alignment', platform='ont',
                 extra_channels=None, custom_depth=None):
        """
        Args:
            mode: 'full_alignment' or 'pileup'
            platform: 'ont', 'hifi', or 'ilmn'
            extra_channels: List of extra channel names to append (e.g., ['dwell_time'])
            custom_depth: Override matrix depth (FA mode only)
        """
        self.mode = mode
        self.platform = platform

        if mode == 'full_alignment':
            base_channels = list(self.FA_CHANNELS)
            self.matrix_depth = custom_depth or self.MATRIX_DEPTH_DICT[platform]
            self.float_type = 'int8'
        elif mode == 'pileup':
            base_channels = list(self.PILEUP_CHANNELS)
            self.matrix_depth = None
            self.float_type = 'int32'
        else:
            raise ValueError("mode must be 'full_alignment' or 'pileup', got: {}".format(mode))

        self.channels = base_channels
        self.extra_channels = []
        if extra_channels:
            self.extra_channels = list(extra_channels)
            self.channels = base_channels + self.extra_channels

        self.no_of_positions = self.NO_OF_POSITIONS
        self.flanking_base_num = self.FLANKING_BASE_NUM
        self.label_size = self.LABEL_SIZE

        # Extra registered components (beyond position_matrix, label, position, alt_info)
        self._extra_components = []

    @property
    def channel_size(self):
        return len(self.channels)

    @property
    def tensor_shape(self):
        if self.mode == 'pileup':
            return [self.no_of_positions, self.channel_size]
        return [self.matrix_depth, self.no_of_positions, self.channel_size]

    def add_component(self, name, shape, dtype='int8', optional=True):
        """
        Register an extra tensor component (stored as additional HDF5 dataset).

        Args:
            name: Dataset name (e.g., 'phasing_matrix', 'signal_features')
            shape: Per-sample shape (e.g., [89, 64])
            dtype: numpy dtype
            optional: Whether this component is optional
        """
        self._extra_components.append(TensorComponent(name, shape, dtype, optional))
        return self

    @property
    def extra_components(self):
        return list(self._extra_components)

    def create_hdf5_datasets(self, h5file, chunk_rows=500):
        """
        Create all HDF5 datasets based on this config.
        Creates the standard 4 datasets + any registered extra components.

        Args:
            h5file: h5py.File object opened in write mode
            chunk_rows: Chunk size for HDF5 compression
        """
        from clair3.utils import _hdf5_compression_kwargs, ensure_hdf5_plugin_path
        ensure_hdf5_plugin_path()
        compression_kwargs = _hdf5_compression_kwargs()

        tensor_shape = tuple(self.tensor_shape)

        # Standard datasets
        h5file.create_dataset(
            "position_matrix",
            shape=(0,) + tensor_shape,
            maxshape=(None,) + tensor_shape,
            chunks=(chunk_rows,) + tensor_shape,
            dtype=np.dtype(self.float_type),
            **compression_kwargs,
        )
        h5file.create_dataset(
            "position",
            shape=(0, 1), maxshape=(None, 1), chunks=(chunk_rows, 1),
            dtype="S{}".format(self.no_of_positions + 50),
            **compression_kwargs,
        )
        h5file.create_dataset(
            "label",
            shape=(0, self.label_size), maxshape=(None, self.label_size),
            chunks=(chunk_rows, self.label_size),
            dtype=np.dtype(self.float_type),
            **compression_kwargs,
        )
        h5file.create_dataset(
            "alt_info",
            shape=(0, 1), maxshape=(None, 1), chunks=(chunk_rows, 1),
            dtype="S5000",
            **compression_kwargs,
        )

        # Extra components
        for comp in self._extra_components:
            comp_shape = tuple(comp.shape)
            h5file.create_dataset(
                comp.name,
                shape=(0,) + comp_shape,
                maxshape=(None,) + comp_shape,
                chunks=(chunk_rows,) + comp_shape,
                dtype=np.dtype(comp.dtype),
                **compression_kwargs,
            )

        return h5file

    def __repr__(self):
        extra = ""
        if self._extra_components:
            extra = ", extra=[{}]".format(", ".join(c.name for c in self._extra_components))
        return "TensorConfig(mode={!r}, platform={!r}, shape={}{})".format(
            self.mode, self.platform, self.tensor_shape, extra)
