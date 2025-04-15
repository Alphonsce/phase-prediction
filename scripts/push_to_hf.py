import os
from datasets import Dataset, Features, Value
from huggingface_hub import login, HfApi

from tqdm import tqdm

TOKEN = "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"

def generate_examples(data_dir):
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"Directory {data_dir} not found")

    for filename in sorted(os.listdir(data_dir)):
        if filename.endswith('.cif'):
            path = os.path.join(data_dir, filename)
            try:
                with open(path, "r", encoding="utf-8") as f:
                    content = f.read()
                yield {
                    "file_name": filename,
                    "content": content
                }
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")
                continue

def create_and_push_dataset(repo_name, data_dir, batch_size=1000):
    login(token=TOKEN)
    features = Features({
        "file_name": Value("string"),
        "content": Value("string"),
    })
    api = HfApi()
    
    batch = []
    shard_num = 0
    for example in tqdm(generate_examples(data_dir), total=len(os.listdir(data_dir)), desc="Uploading to HF..."):
        batch.append(example)
        if len(batch) == batch_size:
            # Create dataset from the current batch
            ds = Dataset.from_list(batch, features=features)
            # Save to a temporary parquet file
            temp_parquet = f"temp_{shard_num}.parquet"
            ds.to_parquet(temp_parquet)
            # Count total number of shards (current + 1 for potential last batch)
            total_shards = (len(os.listdir(data_dir)) // batch_size) + 1
            # Format shard number with leading zeros based on total shards
            formatted_shard = str(shard_num).zfill(5)
            formatted_total = str(total_shards).zfill(5)
            # Upload to the hub with the new naming convention
            path_in_repo = f"data/train-{formatted_shard}-of-{formatted_total}.parquet"
            api.upload_file(
                path_or_fileobj=temp_parquet,
                path_in_repo=path_in_repo,
                repo_id=repo_name,
                repo_type="dataset",
                token=TOKEN,
            )
            # Clean up the temporary file
            os.remove(temp_parquet)
            batch = []
            shard_num += 1
    
    # Process the last batch
    if batch:
        ds = Dataset.from_list(batch, features=features)
        temp_parquet = f"temp_{shard_num}.parquet"
        ds.to_parquet(temp_parquet)
        # Count total number of shards
        total_shards = (len(os.listdir(data_dir)) // batch_size) + 1
        # Format shard number with leading zeros based on total shards
        formatted_shard = str(shard_num).zfill(5)
        formatted_total = str(total_shards).zfill(5)
        # Upload to the hub with the same naming convention as the batches
        path_in_repo = f"data/train-{formatted_shard}-of-{formatted_total}.parquet"
        api.upload_file(
            path_or_fileobj=temp_parquet,
            path_in_repo=path_in_repo,
            repo_id=repo_name,
            repo_type="dataset",
            token=TOKEN,
        )
        os.remove(temp_parquet)

if __name__ == "__main__":
    all_cifs_path = "/Users/aleksandr.varlamov/cif/all_cifs"
    
    create_and_push_dataset("Alphonsce/cif-dataset", all_cifs_path, batch_size=1000)  