
import argparse

def concatBlast(blastItems, seedSize):
    results = {}
    query = blastItems[0]["qseq"]
    for item in blastItems:
        qstart = int(item["qstart"])
        qend = int(item["qend"])
        sstart = int(item["sstart"])
        send = int(item["send"])
        nident = int(item["nident"])
        pident = item["pident"]
        qstrand = "+" if qstart - qend < 0 else "-"
        sstrand = "+" if sstart - send < 0 else "-"
        if qstrand != sstrand:
            # query, subject strand 같아야함
            continue
        if qstrand == "-":
            temp = qstart
            qstart = qend
            qend = temp
        if sstrand == "-":
            temp = sstart
            sstart = send
            send = temp

        if len(results.keys()) == 0:
            results = {
                "qseq": item["qseq"].split(":")[0],
                "sseq": item["sseq"].split(":")[0],
                "qlen": int(item["qlen"]),
                "slen": int(item["slen"]),
                "nident": int(item["nident"]),
                "qstart": qstart,
                "qend": qend,
                "sstart": sstart,
                "send": send,
                "qpos": [str(qstart) + "-" + str(qend)],
                "spos": [str(sstart) + "-" + str(send)],
                "match": [str(nident)],
                "pmatch": [pident],
                "qstrand": qstrand,
                "sstrand": sstrand
            }
        else:
            if qstrand != results["qstrand"] or sstrand != results["sstrand"]:
                continue

            if results["qstart"] > qend:
                if results["send"] <= send:
                    # query left, subject right
                    continue
                if results["sstart"] <= sstart and send <= results["send"]:
                    continue

                if results["sstart"] <= send:
                    nonOverlap = send - results["sstart"] + 1
                    qend -= nonOverlap
                    send -= nonOverlap
                    nident -= nonOverlap

                if nident < seedSize:
                    continue

                results["qstart"] = qstart
                results["sstart"] = sstart
                results["nident"] += nident
                results["qpos"].append(str(qstart) + "-" + str(qend))
                results["spos"].append(str(sstart) + "-" + str(send))
                results["match"].append(str(nident))
                results["pmatch"].append(str(pident))

            elif results["qend"] < qstart:
                if results["sstart"] >= sstart:
                    continue
                if results["sstart"] <= sstart and send <= results["send"]:
                    continue
                if results["send"] >= sstart:
                    nonOverlap = results["send"] - sstart + 1
                    qstart += nonOverlap
                    sstart += nonOverlap
                    nident -= nonOverlap

                if nident < seedSize:
                    continue

                results["qend"] = qend
                results["send"] = send
                results["nident"] += nident
                results["qpos"].append(str(qstart) + "-" + str(qend))
                results["spos"].append(str(sstart) + "-" + str(send))
                results["match"].append(str(nident))
                results["pmatch"].append(str(pident))

            elif results["qstart"] >= qstart and results["qend"] >= qend:
                if results["send"] <= send:
                    continue
                if results["sstart"] <= sstart and send <= results["send"]:
                    continue

                nonOverlap = qend - results["qstart"] + 1
                if results["sstart"] <= send:
                    subjectNontOverlap = send - results["sstart"] + 1
                    nonOverlap = subjectNontOverlap if subjectNontOverlap < nonOverlap else nonOverlap

                if nident < seedSize:
                    continue

                send -= nonOverlap
                qend -= nonOverlap
                nident -= nonOverlap
                results["qstart"] = qstart
                results["sstart"] = sstart
                results["nident"] += nident
                results["qpos"].append(str(qstart) + "-" + str(qend))
                results["spos"].append(str(sstart) + "-" + str(send))
                results["match"].append(str(nident))
                results["pmatch"].append(str(pident))

            elif results["qstart"] < qstart and results["qend"] < qend:
                if results["sstart"] >= sstart:
                    continue
                if results["sstart"] <= sstart and send <= results["send"]:
                    continue

                nonOverlap = results["qend"] - qstart + 1
                if results["send"] >= sstart:
                    subjectNontOverlap = results["send"] - sstart + 1
                    nonOverlap = subjectNontOverlap if subjectNontOverlap < nonOverlap else nonOverlap

                if nident < seedSize:
                    continue

                qstart += nonOverlap
                sstart += nonOverlap
                nident -= nonOverlap
                results["qend"] = qend
                results["send"] = send
                results["nident"] += nident
                results["qpos"].append(str(qstart) + "-" + str(qend))
                results["spos"].append(str(sstart) + "-" + str(send))
                results["match"].append(str(nident))
                results["pmatch"].append(str(pident))

    if "nident" not in results:
        return None
    
    qcov = results["nident"] / results["qlen"]
    scov = results["nident"] / results["slen"]
    normCov = (qcov + scov) / 2
    results["qcov"] = qcov
    results["scov"] = scov
    results["normCov"] = normCov       
   
    return results



def main():
    parser = argparse.ArgumentParser(description="Concatenates fragmented BLAST local alignments into continuous query-subject alignments for accurate coverage and similarity calculations.")

    parser.add_argument(
        "-i", "--input",
        dest="input",
        action="store",
        required=True,
        help="BLAST output file in tabular format (outfmt 6). The input file must be generated using BLAST with '-outfmt 6' to ensure compatibility."
    )

    parser.add_argument(
        "-o", "--output",
        dest="output",
        action="store",
        required=True,
        help="Output file name where the processed results will be saved."
    )

    parser.add_argument(
        "-c", "--coverage",
        dest="coverage",
        action="store",
        default=0.7,
        type=float,
        help="Coverage threshold (default: 0.7). Only alignments with a normalized coverage greater than or equal to this threshold will be included. Range: 0 to 1"
    )

    parser.add_argument(
        "-s", "--seed",
        dest="seed",
        action="store",
        default=23,
        type=int,
        help="Seed size for alignment merging (default: 23). Determines the minimum required overlap to concatenate alignments."
    )

    parser.add_argument(
        "-f", "--format",
        dest="format",
        action="store",
        choices=["result", "link"],
        default="result",
        help=(
            "Specify the output format:\n"
            "  - 'result': Concatenated BLAST results in tab-delimited format (default).\n"
            "  - 'link': JSON-formatted output for visualization, showing linkage information."
        )
    )

    option = parser.parse_args()

    headers = ["qseq", "sseq", "nident", "pident", "positive",
            "ppos", "mismatch", "gaps", "gapopen", "length",
            "qlen", "slen", "qstart", "qend", "sstart",
            "send", "bitscore", "evalue"]

    infile = open(option.input)
    limitCov = option.coverage
    seedSize = option.seed

    blastDict = {}
    queryKeys = {}
    subjectKeys = {}
    for line in infile:
        item = line.strip("\n").split("\t")
        blastDict.setdefault(item[0], []).append({key: value for key, value in zip(headers, item)})
        queryKeys[item[0]] = ""
        subjectKeys[item[1]] = ""
                    
    queryKeys = list(sorted(queryKeys.keys()))
    subjectKeys = list(sorted(subjectKeys.keys()))


    resultDict = {}
    for query in queryKeys:
        blastItems = blastDict[query]
        resultDict[query] = []
        for subject in subjectKeys:
            items = list(filter(lambda x:x["sseq"] == subject, blastItems))
            items = sorted(items, key=lambda x: int(x["nident"]), reverse=True)
            if len(items) == 0:
                continue
            concatResult = concatBlast(items, seedSize)
            if concatResult is not None:
                resultDict[query].append(concatResult)
            


    text = []
    if option.format == "result":
        for query in resultDict:

            items = resultDict[query]
            for item in items:

                if item["qcov"] < limitCov or item["scov"] < limitCov:
                    continue
                result = [item["qseq"], item["sseq"], item["qlen"], item["slen"], item["nident"], item["qstart"],
                        item["qend"], item["sstart"], item["send"], item["qcov"], item["scov"], item["normCov"], ",".join(sorted(item["qpos"], key=lambda x: int(x.split("-")[1]))), ",".join(sorted(item["spos"], key=lambda x: int(x.split("-")[1]))), ",".join(item["match"])]
                text.append("\t".join(map(lambda x: str(x), result)))
    elif option.format == "link":
        for query in resultDict:
            for item in  resultDict[query]:
                if item["qcov"] < limitCov or item["scov"] < limitCov:
                    continue
                for idx in range(len(item["qpos"])):
                    value = list(map(lambda x: str(x), [item["qseq"]] + item["qpos"][idx].split("-") + [item["qcov"]] + [item["sseq"]] + item["spos"][idx].split("-") + [item["scov"], item["match"][idx], item["pmatch"][idx]]))
                    text.append("\t".join(value))

    outfile = open(option.output, "wt")
    outfile.write("\n".join(text))
    outfile.close()



if __name__ == "__main__":
    main()

