//! Fixture física mínima para el lector F4.3, sin habilitar la publicación F4.5.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{build_effective_config, load_and_validate_config};
use genome_population_simulator::biological_contamination::contamination_target_pair_count;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::event_manifest::{
    EventManifest, EventManifestArtifact, EventManifestCoverage, EventManifestDemultiplexing,
    EventManifestReference, event_manifest_bytes, event_manifest_relative_path,
};
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::fastq_validation::F4CoverageReport;
use genome_population_simulator::identity::{
    EventId, IndividualId, LibraryId, ReferenceId, derive_execution_id,
};
use genome_population_simulator::integrity::{configuration_digest, identify_file};
use genome_population_simulator::randomness::{
    F4TechnicalDecision, RootSeed, derive_f4_technical_subseed, f4_uniform_closed,
};
use genome_population_simulator::results::{
    ArtifactRecord, UnitState, validate_published_resources,
};
use genome_population_simulator::technical_identity::{
    TechnicalMate, derive_technical_identity_plan,
};
use serde_json::{Value, json};

static NONCE: AtomicU64 = AtomicU64::new(0);
const ROOT_SEED: &str = "000102030405060708090a0b1c1d1e1f101112131415161718191a1b1c1d1e1f";

#[test]
fn physical_reader_accepts_coherent_f4_documents_and_rejects_independent_chimera_mutations() {
    assert_case("valid", |_| {}, true);
    assert_case(
        "event-rate",
        |fixture| {
            fixture.event["artifacts"]["chimera"]["rate"] = json!(1.0);
        },
        false,
    );
    assert_case(
        "general-rate",
        |fixture| {
            fixture.general["technical_events"][0]["artifacts"]["chimera"]["rate"] = json!(1.0);
        },
        false,
    );
    assert_case(
        "event-counter",
        |fixture| {
            fixture.event["artifacts"]["chimera"]["realized_pairs"] = json!(1);
        },
        false,
    );
    assert_case(
        "general-counter",
        |fixture| {
            fixture.general["technical_events"][0]["artifacts"]["chimera"]["realized_pairs"] =
                json!(1);
        },
        false,
    );
    assert_case(
        "event-index-hopping-rate",
        |fixture| {
            fixture.event["demultiplexing"]["index_hopping_rate"] = json!(0.5);
        },
        false,
    );
    assert_case(
        "event-sequencing-errors-mode",
        |fixture| {
            fixture.event["error_quality"]["requested_mode"] = json!("disabled");
            fixture.event["error_quality"]["effective_mode"] = json!("disabled");
            fixture.event["error_quality"]["qualities_represent_realized_error_probability"] =
                json!(false);
        },
        false,
    );
    for (label, mate, field, value) in [
        ("event-r1-quality-shift", "r1", "quality_shift", json!(1)),
        ("event-r1-error-scale", "r1", "error_scale", json!(2.0)),
        ("event-r2-quality-shift", "r2", "quality_shift", json!(-1)),
        ("event-r2-error-scale", "r2", "error_scale", json!(0.5)),
    ] {
        assert_case(
            label,
            |fixture| {
                fixture.event["error_quality"][mate][field] = value;
            },
            false,
        );
    }
}

#[test]
fn physical_reader_rederives_the_neutral_contamination_summary_and_rejects_mutations() {
    assert_contamination_case("valid", |_| {}, true);
    assert_contamination_case(
        "recipient-library",
        |fixture| {
            fixture.event["artifacts"]["biological_contamination"][0]["recipient_library_id"] =
                json!("LIB_wrong");
        },
        false,
    );
    assert_contamination_case(
        "range",
        |fixture| {
            fixture.event["artifacts"]["biological_contamination"][0]["percent_range"]["maximum"] =
                json!(0.5);
        },
        false,
    );
    assert_contamination_case(
        "rounded-target",
        |fixture| {
            fixture.event["artifacts"]["biological_contamination"][0]["target_pair_count"] =
                json!(1);
        },
        false,
    );
    assert_contamination_case(
        "realized-counts",
        |fixture| {
            fixture.event["artifacts"]["biological_contamination"][0]["retained_own_pair_count"] =
                json!(4);
        },
        false,
    );
    assert_contamination_case(
        "percent-subseed",
        |fixture| {
            fixture.event["artifacts"]["biological_contamination"][0]["subseeds"]["percent"] =
                json!("f".repeat(64));
        },
        false,
    );
    assert_contamination_sampled_percent_rng_guard();
}

fn assert_case(label: &str, mutate: impl FnOnce(&mut Fixture), valid: bool) {
    let mut fixture = Fixture::new(label);
    fixture.event["artifacts"]["biological_contamination"] =
        json!([neutral_contamination_summary(&fixture)]);
    mutate(&mut fixture);
    fixture.persist();
    let result = validate_published_resources(
        &fixture.output,
        &fixture.effective_path,
        &fs::read(&fixture.effective_path).expect("effective bytes"),
        &fixture.manifest_path,
        &fs::read(&fixture.manifest_path).expect("general bytes"),
    );
    if valid {
        result.expect("fixture física F4 coherente");
    } else {
        assert!(result.is_err(), "{label} debe ser rechazado");
    }
    fs::remove_dir_all(fixture.root).expect("cleanup");
}

fn assert_contamination_case(label: &str, mutate: impl FnOnce(&mut Fixture), valid: bool) {
    let mut fixture = Fixture::new(label);
    fixture.event["artifacts"]["biological_contamination"] =
        json!([neutral_contamination_summary(&fixture)]);
    mutate(&mut fixture);
    fixture.persist();
    let result = validate_published_resources(
        &fixture.output,
        &fixture.effective_path,
        &fs::read(&fixture.effective_path).expect("effective bytes"),
        &fixture.manifest_path,
        &fs::read(&fixture.manifest_path).expect("general bytes"),
    );
    if valid {
        result.expect("resumen neutro coherente");
    } else {
        assert!(result.is_err(), "{label} debe ser rechazado");
    }
    fs::remove_dir_all(fixture.root).expect("cleanup");
}

fn assert_contamination_sampled_percent_rng_guard() {
    let mut fixture = Fixture::new("sampled-percent-rng-guard");
    fixture.event["artifacts"]["biological_contamination"] =
        json!([neutral_contamination_summary(&fixture)]);
    let record = &mut fixture.event["artifacts"]["biological_contamination"][0];
    let sampled = record["sampled_percent"]
        .as_f64()
        .expect("sampled percentage");
    let mutated = if sampled.to_bits() == 0.0_f64.to_bits() {
        1.0
    } else {
        0.0
    };
    let clean = record["clean_pair_count"]
        .as_u64()
        .expect("clean pair count");
    let target = contamination_target_pair_count(clean, mutated).expect("coherent target");
    assert_eq!(target, 0, "the low non-degenerate range needs no donors");
    record["sampled_percent"] = json!(mutated);
    record["target_pair_count"] = json!(target);
    record["donor_contributions"] = json!([]);
    fixture.persist();

    let error = validate_published_resources(
        &fixture.output,
        &fixture.effective_path,
        &fs::read(&fixture.effective_path).expect("effective bytes"),
        &fixture.manifest_path,
        &fs::read(&fixture.manifest_path).expect("general bytes"),
    )
    .expect_err("a different in-range percentage must fail the RNG bit guard");
    assert!(
        error
            .render_human()
            .contains("porcentaje contaminante no coincide con la subsemilla F4"),
        "the rejection must come from the sampled-percent RNG/to_bits guard: {}",
        error.render_human()
    );
    fs::remove_dir_all(fixture.root).expect("cleanup");
}

fn neutral_contamination_summary(fixture: &Fixture) -> Value {
    let recipient = fixture.event["selected_individuals"][0]["individual_id"]
        .as_str()
        .expect("recipient");
    let library = fixture.event["selected_individuals"][0]["library_id"]
        .as_str()
        .expect("library");
    let root = RootSeed::parse("fixture root", ROOT_SEED).expect("root seed");
    let event = EventId::parse("reads-main", "fixture").expect("event");
    let recipient = IndividualId::parse(recipient, "fixture").expect("recipient ID");
    let library = LibraryId::parse(library, "fixture").expect("library ID");
    let percent_subseed = derive_f4_technical_subseed(
        root,
        F4TechnicalDecision::contamination_percent(&event, &library, &recipient),
    )
    .expect("percent subseed");
    let sampled_percent =
        f4_uniform_closed(&mut percent_subseed.local_rng(), 0.0, 1.0).expect("sampled percent");
    let target_pair_count =
        contamination_target_pair_count(5, sampled_percent).expect("contamination target");
    assert_eq!(target_pair_count, 0, "the fixture has no donor library");
    let seed = |decision| {
        derive_f4_technical_subseed(root, decision)
            .expect("subseed")
            .to_hex()
    };
    json!({
        "recipient_id": recipient.as_str(),
        "recipient_library_id": library.as_str(),
        "percent_range": {"minimum": 0.0, "maximum": 1.0},
        "sampled_percent": sampled_percent,
        "clean_pair_count": 5,
        "target_pair_count": target_pair_count,
        "removed_own_pair_count": 0,
        "inserted_contaminant_pair_count": 0,
        "retained_own_pair_count": 5,
        "final_pair_count": 5,
        "subseeds": {
            "percent": percent_subseed.to_hex(),
            "donors": seed(F4TechnicalDecision::contamination_donors(&event, &library, &recipient)),
            "removal": seed(F4TechnicalDecision::contamination_recipients(&event, &library, &recipient)),
        },
        "donor_contributions": [],
    })
}

struct Fixture {
    root: PathBuf,
    output: PathBuf,
    effective_path: PathBuf,
    manifest_path: PathBuf,
    event_path: PathBuf,
    event_artifact_id: String,
    event: Value,
    general: Value,
}

impl Fixture {
    fn new(label: &str) -> Self {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("fixture root");
        let source = root.join("general.toml");
        fs::write(&source, base_toml()).expect("base configuration");
        let output = root.join("output");
        execute_validated_config(
            &load_and_validate_config(&source).expect("baseline configuration"),
            &output,
        )
        .expect("baseline publication without F4 event");

        let technical_toml = technical_toml();
        let effective_path = output.join("effective-config.json");
        let manifest_path = output.join("manifest.json");
        let mut general = read_json(&manifest_path);
        let mut effective = read_json(&effective_path);
        let technical = technical_effective_event(&technical_toml);
        effective["technical_event"] = technical.clone();
        write_json(&effective_path, &effective);
        refresh_catalog_record(&output, &mut general, "effective-config");

        let configuration_id = general["configuration_input"]["artifact_id"]
            .as_str()
            .expect("configuration artifact id")
            .to_owned();
        let configuration_path = output.join(catalog_path(&general, &configuration_id));
        fs::write(&configuration_path, &technical_toml).expect("physical configuration");
        refresh_catalog_record(&output, &mut general, &configuration_id);
        let root_seed = RootSeed::parse("fixture root", ROOT_SEED).expect("root seed");
        let digest = configuration_digest(technical_toml.as_bytes()).expect("configuration digest");
        general["execution_id"] = json!(
            derive_execution_id(root_seed, &digest)
                .expect("execution id")
                .as_str()
        );

        let reference = artifact_record(
            general["artifacts"]
                .as_array()
                .expect("catalog")
                .iter()
                .find(|artifact| artifact["kind"] == "reference_fasta")
                .expect("reference artifact"),
        );
        let event_config =
            build_effective_config(&declared(&technical_toml)).expect("event config");
        let event = event_config.technical_event().expect("technical event");
        let reference_id = ReferenceId::parse("reference-main", "fixture").expect("reference id");
        let identities = derive_technical_identity_plan(
            root_seed,
            &reference_id,
            event.event_id().clone(),
            event.individual_ids(),
        )
        .expect("technical identities");
        let library = identities
            .library(event.individual_ids().first().expect("selected individual"))
            .expect("library identity");
        let (r1_relative, r2_relative) =
            genome_population_simulator::fastq_output::paired_fastq_relative_paths(
                event.event_id().as_str(),
                library.individual_id().as_str(),
            );
        let r1 = write_artifact(&output, "fastq-LIB_sample-A", "fastq", &r1_relative);
        let r2 = write_artifact(&output, "fastq-LIB_sample-A-r2", "fastq", &r2_relative);
        let outputs = [
            EventManifestArtifact {
                id: &r1.id,
                kind: "fastq",
                state: UnitState::Complete,
                individual_id: library.individual_id(),
                library_id: library.library_id(),
                mate: TechnicalMate::R1,
                artifact: Some(&r1),
                diagnostic: None,
            },
            EventManifestArtifact {
                id: &r2.id,
                kind: "fastq",
                state: UnitState::Complete,
                individual_id: library.individual_id(),
                library_id: library.library_id(),
                mate: TechnicalMate::R2,
                artifact: Some(&r2),
                diagnostic: None,
            },
        ];
        let coverage = [EventManifestCoverage {
            individual_id: library.individual_id(),
            library_id: library.library_id(),
            report: F4CoverageReport::from_base_counts(1.0, 300, 300, 0, 100, 0).expect("coverage"),
        }];
        let summary = [
            genome_population_simulator::event_manifest::EventManifestFinalFastqSummary {
                individual_id: library.individual_id(),
                library_id: library.library_id(),
                pair_count: 1,
                read_count: 2,
                base_count: 300,
                own_pair_count: 1,
                contaminant_pair_count: 0,
                own_base_count: 300,
                contaminant_base_count: 0,
                n_base_count: 0,
            },
        ];
        let document = EventManifest {
            event,
            state: UnitState::Complete,
            identities: &identities,
            reference: EventManifestReference {
                reference_id: &reference_id,
                fasta: &reference,
            },
            artifacts: &outputs,
            coverage: &coverage,
            final_fastq_summary: &summary,
            error_quality:
                genome_population_simulator::event_manifest::EventManifestErrorQuality::neutral(),
            demultiplexing: EventManifestDemultiplexing::disabled(),
            chimera_rate: event.chimera_rate(),
            realized_chimera_pairs: 0,
            biological_contamination: &[],
            stutter: None,
            diagnostics: &[],
        };
        let event_path = output.join(event_manifest_relative_path(event.event_id()));
        fs::create_dir_all(event_path.parent().expect("event directory")).expect("event directory");
        fs::write(
            &event_path,
            event_manifest_bytes(&document).expect("event bytes"),
        )
        .expect("event manifest");
        let event_record = ArtifactRecord::from_identity(
            "event-reads-main-manifest",
            "event_manifest",
            &identify_file(&event_path, Path::new("events/reads-main/manifest.json"))
                .expect("event identity"),
        )
        .expect("event record");
        append_catalog(
            &mut general,
            &[r1.clone(), r2.clone(), event_record.clone()],
        );
        general["events"] = json!([{
            "event_id": event.event_id().as_str(),
            "state": "complete",
            "manifest_artifact_id": event_record.id,
        }]);
        let effective_record = catalog_record(&general, "effective-config");
        general["technical_events"] = json!([{
            "event_id": event.event_id().as_str(),
            "effective_config": effective_record,
            "profile": {"id": "gps-hiseqx-pcrfree-pe150-v1", "sha256": "3d8878ab920237dae03705276a26dbb0823941830d531374d35aa1c540f6deea"},
            "rng_catalog": "f4-technical-v3",
            "id_catalog": {"version":"3", "derivations":{"v1":["execution","reference","individual","event","batch","library","read_group","template"],"v2":["PreHopInstanceId","TechnicalPairId"]}},
            "artifacts": {"chimera": {"rate": 0.0, "realized_pairs": 0}},
            "event_manifest": event_record.to_json(),
        }]);
        general["output_categories"]["minimum"] = general["artifacts"]
            .as_array()
            .expect("catalog")
            .iter()
            .map(|artifact| artifact["id"].clone())
            .collect();

        Self {
            root,
            output,
            effective_path,
            manifest_path,
            event_path: event_path.clone(),
            event_artifact_id: event_record.id,
            event: serde_json::from_slice(&fs::read(&event_path).expect("event bytes"))
                .expect("event JSON"),
            general,
        }
    }

    fn persist(&mut self) {
        write_json(&self.event_path, &self.event);
        refresh_catalog_record(&self.output, &mut self.general, &self.event_artifact_id);
        self.general["technical_events"][0]["event_manifest"] =
            catalog_record(&self.general, &self.event_artifact_id);
        write_json(&self.manifest_path, &self.general);
    }
}

fn append_catalog(general: &mut Value, records: &[ArtifactRecord]) {
    let catalog = general["artifacts"].as_array_mut().expect("catalog");
    catalog.extend(records.iter().map(ArtifactRecord::to_json));
}

fn write_artifact(root: &Path, id: &str, kind: &str, relative: &Path) -> ArtifactRecord {
    let path = root.join(relative);
    fs::create_dir_all(path.parent().expect("artifact directory")).expect("artifact directory");
    fs::write(&path, b"@fixture\nA\n+\nI\n").expect("artifact bytes");
    ArtifactRecord::from_identity(id, kind, &identify_file(&path, relative).expect("identity"))
        .expect("artifact record")
}

fn refresh_catalog_record(root: &Path, general: &mut Value, id: &str) {
    let relative = catalog_path(general, id);
    let identity = identify_file(&root.join(&relative), Path::new(&relative)).expect("identity");
    let record = general["artifacts"]
        .as_array_mut()
        .expect("catalog")
        .iter_mut()
        .find(|artifact| artifact["id"] == id)
        .expect("catalog record");
    record["size_bytes"] = json!(identity.size_bytes());
    record["checksum"]["digest"] = json!(identity.digest());
}

fn catalog_path(general: &Value, id: &str) -> String {
    general["artifacts"]
        .as_array()
        .expect("catalog")
        .iter()
        .find(|artifact| artifact["id"] == id)
        .and_then(|artifact| artifact["path"].as_str())
        .expect("catalog path")
        .to_owned()
}

fn catalog_record(general: &Value, id: &str) -> Value {
    general["artifacts"]
        .as_array()
        .expect("catalog")
        .iter()
        .find(|artifact| artifact["id"] == id)
        .expect("catalog record")
        .clone()
}

fn artifact_record(value: &Value) -> ArtifactRecord {
    ArtifactRecord {
        id: value["id"].as_str().expect("id").to_owned(),
        kind: value["kind"].as_str().expect("kind").to_owned(),
        path: value["path"].as_str().expect("path").to_owned(),
        size_bytes: value["size_bytes"].as_u64().expect("size"),
        digest: value["checksum"]["digest"]
            .as_str()
            .expect("digest")
            .to_owned(),
    }
}

fn technical_effective_event(source: &str) -> Value {
    let effective = build_effective_config(&declared(source)).expect("effective configuration");
    let reference_id = ReferenceId::parse("reference-main", "fixture").expect("reference id");
    let reference = ArtifactRecord {
        id: "reference-main-fasta".to_owned(),
        kind: "reference_fasta".to_owned(),
        path: "references/reference-main.fasta".to_owned(),
        size_bytes: 1,
        digest: "a".repeat(64),
    };
    let mut reference_model = genome_population_simulator::reference::Reference::new(
        reference_id,
        genome_population_simulator::reference::ReferenceOrigin::Generated,
        vec![
            genome_population_simulator::reference::ReferenceSequence::new(
                genome_population_simulator::identity::SequenceId::parse("chr1", "fixture")
                    .expect("sequence id"),
                genome_population_simulator::reference::FastaHeader::parse(b"chr1")
                    .expect("header"),
                genome_population_simulator::reference::SequenceClass::Chromosome,
                &genome_population_simulator::reference::DnaSequence::try_from_bytes(&[b'A'; 100])
                    .expect("DNA"),
            )
            .expect("sequence"),
        ],
    )
    .expect("reference");
    reference_model
        .bind_authoritative_resource(reference.path.clone().into())
        .expect("resource");
    genome_population_simulator::results::effective_config_document_with_reference(
        &effective,
        genome_population_simulator::results::ReferenceOutput {
            reference: &reference_model,
            fasta: &reference,
        },
    )
    .expect("effective document")["technical_event"]
        .clone()
}

fn declared(source: &str) -> genome_population_simulator::configuration::DeclaredGeneralConfig {
    parse_declared_config(
        admit_toml_header(
            Path::new("f4-physical-fixture.toml"),
            source.as_bytes().to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared configuration")
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("JSON bytes")).expect("JSON")
}

fn write_json(path: &Path, document: &Value) {
    fs::write(path, serde_json::to_vec(document).expect("JSON")).expect("write JSON");
}

fn base_toml() -> String {
    format!(
        "config_kind = \"general\"\nschema_version = \"2\"\nseed = \"{ROOT_SEED}\"\n\n[reference]\nkind = \"simulated\"\nid = \"reference-main\"\n\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\n[[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 100\n\n[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 2\ngenotype_model = \"hwe\"\nindividual_ids = [\"sample-A\"]\n\n[population_snvs]\nmode = \"generate\"\n\n[population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\nspatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n\n[[population_snvs.generate.populations]]\npopulation_id = \"P1\"\nfrequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n"
    )
}

fn technical_toml() -> String {
    format!(
        "{}\n[technical_event]\nevent_id = \"reads-main\"\nindividual_ids = [\"sample-A\"]\ntarget_coverage_wgs = 1\ncontamination_min_percent = 0\ncontamination_max_percent = 1\n\n[[technical_event.libraries]]\nindividual_id = \"sample-A\"\ni7 = \"ACGTACGT\"\ni5 = \"TGCATGCA\"\n",
        base_toml()
    )
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-physical-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
