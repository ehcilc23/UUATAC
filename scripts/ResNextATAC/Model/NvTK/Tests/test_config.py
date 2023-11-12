import json
from NvTK.config import *

config = load_config_from_json("./config.json")

model = get_model_from_config(config)

optimizer = get_optimizer_from_config(config, model)

criterion = get_criterion_from_config(config)

trainer_args = parse_trainer_args(config)
