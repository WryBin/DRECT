import sys
import argparse
import pytorch_lightning as pl

from pytorch_lightning.callbacks import ModelCheckpoint


def alter(file,old_str,new_str):
    """
    Replace characters in a file
    """
    file_data = ""
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if old_str in line:
                line = line.replace(old_str,new_str)
            file_data += line
    with open(file,"w",encoding="utf-8") as f:
        f.write(file_data)
 

if __name__ == '__main__':

    pl.seed_everything(42)

    parser = argparse.ArgumentParser()
    parser.add_argument('--Type', type=str, default='DOSY')
    args = parser.parse_args()

    # Change to the correct type
    if args.Type not in ['DOSY', 'T1T2', 'VD']:
        print("Get wrong Type, please check again ('DOSY', 'T1T2', 'VD')")
        sys.exit(1)
    print("Type: {}".format(args.Type))
    
    alter("config.py", "='T1T2'", "='"+args.Type+"'")
    alter("config.py", "='VD'", "='"+args.Type+"'")
    alter("config.py", "='DOSY'", "='"+args.Type+"'")

    import config
    import dataset
    from model import make_model

    # data
    train_loader, val_loader = dataset.load_dataloader(batch_size=config.batch_size)

    # initial model
    model = make_model()

    checkpoint_callback = ModelCheckpoint(dirpath='Result/' + config.Type + '/', monitor='val_loss', save_last=True, save_top_k=100)

    trainer = pl.Trainer(accelerator='gpu', gpus=[int(config.gpu)], max_epochs=config.max_epochs, callbacks=checkpoint_callback)

    trainer.fit(model, train_loader, val_loader)