use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{anyhow, Context, Result};
use clap::arg_enum;
use log::info;
use rayon::prelude::*;
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{
    max_hash_for_scaled, HashFunctions, KmerMinHash, KmerMinHashBTree,
};
use sourmash::sketch::Sketch;
use structopt::StructOpt;

arg_enum! {
    #[derive(Debug)]
    enum Encodings {
        Protein,
        Hp,
        Dayhoff,
    }
}

impl From<Encodings> for HashFunctions {
    fn from(e: Encodings) -> HashFunctions {
        match e {
            Encodings::Protein => HashFunctions::murmur64_protein,
            Encodings::Hp => HashFunctions::murmur64_hp,
            Encodings::Dayhoff => HashFunctions::murmur64_dayhoff,
        }
    }
}

#[derive(StructOpt, Debug)]
struct Cli {
    /// Query to be subtracted
    #[structopt(parse(from_os_str))]
    query: PathBuf,

    /// List of signatures to remove the query from
    #[structopt(parse(from_os_str))]
    siglist: PathBuf,

    /// ksize
    #[structopt(short = "k", long = "ksize", default_value = "31")]
    ksize: u8,

    /// scaled
    #[structopt(short = "s", long = "scaled", default_value = "10")]
    scaled: usize,

    /// sequence encoding type
    #[structopt(short = "e", long = "encoding", possible_values = &Encodings::variants(), case_insensitive = true, default_value = "Protein")]
    encoding: Encodings,

    /// The path for output
    #[structopt(parse(from_os_str), short = "o", long = "output")]
    output: Option<PathBuf>,
}

fn subtract<P: AsRef<Path>>(
    query: P,
    siglist: P,
    ksize: u8,
    scaled: usize,
    hash_function: HashFunctions,
    output: Option<P>,
) -> Result<()> {
    info!("Loading queries");

    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .hash_function(hash_function)
        .max_hash(max_hash)
        .build();
    let template = Sketch::MinHash(template_mh);

    let query_sig = Signature::from_path(query).unwrap();
    let mut query = None;
    for sig in &query_sig {
        if let Some(sketch) = sig.select_sketch(&template) {
            if let Sketch::MinHash(mh) = sketch {
                query = Some(mh.clone());
            }
        }
    }
    let query = query.unwrap_or_else(|| {
        panic!(
            "Unable to load a sketch matching the provided template: {:?}",
            &template
        )
    });

    info!("Loaded query signature, k={}", ksize);
    let hashes_to_remove = query.mins();

    info!("Loading siglist");
    let siglist_file = BufReader::new(File::open(siglist)?);
    let search_sigs: Vec<PathBuf> = siglist_file
        .lines()
        .map(|line| {
            let mut path = PathBuf::new();
            path.push(line.unwrap());
            path
        })
        .collect();
    info!("Loaded {} sig paths in siglist", search_sigs.len());

    let mut outdir: PathBuf = if let Some(p) = output {
        p.as_ref().into()
    } else {
        let mut path = PathBuf::new();
        path.push("outputs");
        path
    };
    outdir.push(format!("{}", ksize));
    std::fs::create_dir_all(&outdir)?;

    let processed_sigs = AtomicUsize::new(0);

    search_sigs.par_iter().try_for_each(|filename| {
        let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
        if i % 1000 == 0 {
            info!("Processed {} sigs", i);
        }

        let mut search_sig = Signature::from_path(&filename)?.swap_remove(0);

        let mut search_mh = select_and_downsample(&search_sig, &template)
            .with_context(|| format!("Unable to load a sketch from {:?}", &filename,))?;
        // remove the hashes
        search_mh.remove_many(&hashes_to_remove)?;

        // save to output dir
        let mut path = outdir.clone();
        path.push(
            filename
                .file_name()
                .ok_or_else(|| anyhow!("Error converting filename"))?,
        );

        let mut out = BufWriter::new(File::create(path)?);
        search_sig.reset_sketches();
        search_sig.push(Sketch::LargeMinHash(search_mh));
        serde_json::to_writer(&mut out, &[search_sig])?;

        Ok(())
    })
}

fn select_and_downsample(search_sig: &Signature, template: &Sketch) -> Result<KmerMinHashBTree> {
    let mut search_mh: Option<KmerMinHashBTree> = None;
    if let Sketch::MinHash(template) = template {
        for sketch in search_sig.sketches() {
            if let Sketch::MinHash(mh) = sketch {
                if mh.check_compatible(&template).is_ok() {
                    search_mh = Some(mh.into());
                } else if mh.ksize() == template.ksize()
                    && mh.hash_function() == template.hash_function()
                    && mh.seed() == template.seed()
                    && mh.scaled() < template.scaled()
                {
                    search_mh = Some(mh.downsample_max_hash(template.max_hash())?.into());
                }
            }
        }
    }

    search_mh.ok_or_else(|| anyhow!("No sketch matching the provided template: {:?}", &template))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let opts = Cli::from_args();

    subtract(
        opts.query,
        opts.siglist,
        opts.ksize,
        opts.scaled,
        opts.encoding.into(),
        opts.output,
    )?;

    Ok(())
}
