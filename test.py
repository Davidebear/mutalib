import torch
# Copy and paste below code to test CUDA install and determine device
print(f"CUDA present?:                   {torch.cuda.is_available()}")
print(f"CUDA version:                    {torch.version.cuda}")
  
cuda_id = torch.cuda.current_device()
print(f"Identity of current CUDA device: {torch.cuda.current_device()}")
        
print(f"Name of CUDA device:             {torch.cuda.get_device_name(cuda_id)}")