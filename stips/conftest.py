# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure.

import os
import pytest


if not os.environ.get('DATA_TEST', False):
    @pytest.fixture(scope='session')
    def data_base():
        from stips import stips_data_base
        return os.path.join(stips_data_base, 'test')
