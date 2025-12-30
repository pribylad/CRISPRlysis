# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name
# %%
import time
from datetime import datetime

# Get current date and time
now = datetime.now()
# Format the date to YYMMDD
date = now.strftime('%y%m%d')
start_time = time.time()
now
# %%
end_time = time.time()
total_time = end_time - start_time
total_time