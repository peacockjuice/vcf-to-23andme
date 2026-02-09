"""Built-in database of known genetic markers."""

from __future__ import annotations

KNOWN_MARKERS: dict[str, dict[str, str]] = {
    "rs429358": {
        "gene": "APOE",
        "description": "Associated with Alzheimer's disease risk",
    },
    "rs7412": {
        "gene": "APOE",
        "description": "APOE variant affecting lipid profile",
    },
    "rs1801133": {
        "gene": "MTHFR",
        "description": "Associated with folate metabolism",
    },
    "rs1799971": {
        "gene": "OPRM1",
        "description": "May affect pain sensitivity",
    },
    "rs671": {
        "gene": "ALDH2",
        "description": "Associated with alcohol metabolism",
    },
    "rs6265": {
        "gene": "BDNF",
        "description": "Affects neurotrophic activity",
    },
    "rs4680": {
        "gene": "COMT",
        "description": "Involved in dopamine metabolism",
    },
    "rs9939609": {
        "gene": "FTO",
        "description": "Associated with obesity predisposition",
    },
    "rs1805007": {
        "gene": "MC1R",
        "description": "Associated with pigmentation and hair color",
    },
    "rs12913832": {
        "gene": "HERC2",
        "description": "Associated with eye color",
    },
    "rs1229984": {
        "gene": "ADH1B",
        "description": "Affects alcohol metabolism",
    },
    "rs2228145": {
        "gene": "IL6R",
        "description": "Associated with inflammatory processes",
    },
    "rs7903146": {
        "gene": "TCF7L2",
        "description": "Associated with type 2 diabetes risk",
    },
    "rs1801131": {
        "gene": "MTHFR",
        "description": "Alternative variant associated with folate metabolism",
    },
    "rs662": {
        "gene": "PON1",
        "description": "Involved in antioxidant defense",
    },
    "rs1042713": {
        "gene": "ADRB2",
        "description": "Associated with stress response and receptor function",
    },
    "rs1800562": {
        "gene": "HFE",
        "description": "Associated with hemochromatosis",
    },
    "rs1800795": {
        "gene": "IL6",
        "description": "Associated with inflammatory processes",
    },
}
