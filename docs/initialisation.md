[← Back to main page](index.md)

# Initialisation

After successful installation of SeqNado ([Installation](installation.md)).

## Initialize the Workflow
Run the following command to initialize the SeqNado workflow in your project directory:

```bash
seqnado init
```

This will create the necessary configuration files and directory structure for your environment.

## Additional Initialization Options

For all available flags and details, see the CLI reference: [seqnado init](cli.md#cli-seqnado-init).

### Example Usage

To initialize SeqNado with preset genomes and verbose logging:

```bash
seqnado init --preset --verbose
```

To perform a dry run and preview the initialization steps:

```bash
seqnado init --dry-run
```