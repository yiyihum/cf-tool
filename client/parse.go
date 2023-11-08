package client

import (
	"bytes"
	"encoding/json"
	"fmt"
	"html"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"sync"

	"github.com/yiyihum/cf-tool/util"

	"github.com/k0kubun/go-ansi"

	"github.com/fatih/color"
)

func findSample(body []byte) (input [][]byte, output [][]byte, err error) {
	irg := regexp.MustCompile(`class="input"[\s\S]*?<pre>([\s\S]*?)</pre>`)
	org := regexp.MustCompile(`class="output"[\s\S]*?<pre>([\s\S]*?)</pre>`)
	a := irg.FindAllSubmatch(body, -1)
	b := org.FindAllSubmatch(body, -1)
	if a == nil || b == nil || len(a) != len(b) {
		return nil, nil, fmt.Errorf("cannot parse sample with input %v and output %v", len(a), len(b))
	}
	// adapt new codeforces input format
	regexp1, _ := regexp.Compile(`<div(?U).*>`)
	regexp2, _ := regexp.Compile(`</div>`)
	for i := 0; i < len(a); i++ {
		str := regexp1.ReplaceAll(a[i][1], []byte(""))
		a[i][1] = regexp2.ReplaceAll(str, []byte("\n"))
	}
	newline := regexp.MustCompile(`<[\s/br]+?>`)
	filter := func(src []byte) []byte {
		src = newline.ReplaceAll(src, []byte("\n"))
		s := html.UnescapeString(string(src))
		return []byte(strings.TrimSpace(s) + "\n")
	}
	for i := 0; i < len(a); i++ {
		input = append(input, filter(a[i][1]))
		output = append(output, filter(b[i][1]))
	}
	return
}

func findStatement(body []byte) (statementJSON []byte, err error) {
	re_time, _ := regexp.Compile(`time limit per test</div>(.*?)</div>`)
	time := re_time.FindSubmatch(body)
	re_memory, _ := regexp.Compile(`memory limit per test</div>(.*?)</div>`)
	memory := re_memory.FindSubmatch(body)
	standardIO := true
	if !bytes.Contains(body, []byte(`<div class="input-file"><div class="property-title">input</div>standard input</div><div class="output-file"><div class="property-title">output</div>standard output</div>`)) {
		standardIO = false
	}
	rg_problem, _ := regexp.Compile(`</div></div><div><p>([\s\S]*?)</div><div class="input-specification">`)
	problem := rg_problem.FindSubmatch(body)
	rg_input, _ := regexp.Compile(`Input</div><p>([\s\S]*?)</div><div class="output-specification">`)
	inputFormat := rg_input.FindSubmatch(body)
	rg_output, _ := regexp.Compile(`Output</div><p>([\s\S]*?)</div><div class="sample-tests">`)
	outputFormat := rg_output.FindSubmatch(body)
	rg_note, _ := regexp.Compile(`Note</div><p>([\s\S]*?)</p></div></div><p>`)
	note := rg_note.FindSubmatch(body)
	if len(note) < 2 {
		note = [][]byte{[]byte("none"), []byte("none")}
	}
	//rg_tags, _ := regexp.Compile(`<span class="tag-box" style="font-size:1.2rem;" title="Sortings, orderings">\s*(.*?)\s*</span>`)
	//tags := rg_tags.FindAllSubmatch(body, -1)
	//tag := ""
	//if len(tags) == 0 {
	//	tag = "none"
	//} else {
	//	tag = string(tags[0][1])
	//	for i := 1; i < len(tags); i++ {
	//		tag = tag + ", " + string(tags[i][1])
	//	}
	//}
	input, output, err := findSample(body)
	if err != nil {
		return
	}
	examples := make([]map[string]string, len(input))
	for i := 0; i < len(input); i++ {
		examples[i] = map[string]string{
			"input":  string(input[i]),
			"output": string(output[i]),
		}
	}

	newline := regexp.MustCompile(`</?(br|p|div|li|ul)[^>]*>`)
	filter := func(src []byte) string {
		src = newline.ReplaceAll(src, []byte("\n"))
		s := html.UnescapeString(string(src))
		s = strings.TrimSpace(s)
		re := regexp.MustCompile(`\n+`)
		s = re.ReplaceAllString(s, "\n")
		return s
	}

	statementMap := map[string]interface{}{
		"time limit":    filter(time[1]),
		"memory limit":  filter(memory[1]),
		"standard IO":   standardIO,
		"problem":       filter(problem[1]),
		"input format":  filter(inputFormat[1]),
		"output format": filter(outputFormat[1]),
		"note":          filter(note[1]),
		"examples":      examples,
	}
	statementJSON, err = json.Marshal(statementMap)
	if err != nil {
		return nil, err
	}
	return statementJSON, nil
}

// ParseProblem parse problem to path. mu can be nil
func (c *Client) ParseProblem(URL, path string, mu *sync.Mutex) (samples int, standardIO bool, err error) {
	body, err := util.GetBody(c.client, URL)
	if err != nil {
		return
	}

	_, err = findHandle(body)
	if err != nil {
		return
	}
	// save problem body
	fileBody := filepath.Join(path, "problem_body.html")
	e1 := os.WriteFile(fileBody, body, 0644)
	if e1 != nil {
		if mu != nil {
			mu.Lock()
		}
		color.Red(e1.Error())
		if mu != nil {
			mu.Unlock()
		}
	}

	input, output, err := findSample(body)
	if err != nil {
		return
	}

	statementJSON, err := findStatement(body)
	if err != nil {
		return
	}

	standardIO = true
	if !bytes.Contains(body, []byte(`<div class="input-file"><div class="property-title">input</div>standard input</div><div class="output-file"><div class="property-title">output</div>standard output</div>`)) {
		standardIO = false
	}

	for i := 0; i < len(input); i++ {
		fileIn := filepath.Join(path, fmt.Sprintf("in%v.txt", i+1))
		fileOut := filepath.Join(path, fmt.Sprintf("ans%v.txt", i+1))
		e := os.WriteFile(fileIn, input[i], 0644)
		if e != nil {
			if mu != nil {
				mu.Lock()
			}
			color.Red(e.Error())
			if mu != nil {
				mu.Unlock()
			}
		}
		e = os.WriteFile(fileOut, output[i], 0644)
		if e != nil {
			if mu != nil {
				mu.Lock()
			}
			color.Red(e.Error())
			if mu != nil {
				mu.Unlock()
			}
		}
	}
	fileStatement := filepath.Join(path, "statement.json")
	e := os.WriteFile(fileStatement, statementJSON, 0644)
	if e != nil {
		if mu != nil {
			mu.Lock()
		}
		color.Red(e.Error())
		if mu != nil {
			mu.Unlock()
		}
	}
	return len(input), standardIO, nil
}

// Parse parse
func (c *Client) Parse(info Info) (problems []string, paths []string, err error) {
	color.Cyan("Parse " + info.Hint())

	problemID := info.ProblemID
	info.ProblemID = "%v"
	urlFormatter, err := info.ProblemURL(c.host)
	if err != nil {
		return
	}
	info.ProblemID = ""
	if problemID == "" {
		statics, err := c.Statis(info)
		if err != nil {
			return nil, nil, err
		}
		problems = make([]string, len(statics))
		for i, problem := range statics {
			problems[i] = problem.ID
		}
	} else {
		problems = []string{problemID}
	}
	contestPath := info.Path()
	ansi.Printf(color.CyanString("The problem(s) will be saved to %v\n"), color.GreenString(contestPath))

	wg := sync.WaitGroup{}
	wg.Add(len(problems))
	mu := sync.Mutex{}
	paths = make([]string, len(problems))
	for i, problemID := range problems {
		paths[i] = filepath.Join(contestPath, strings.ToLower(problemID))
		go func(problemID, path string) {
			defer wg.Done()
			mu.Lock()
			fmt.Printf("Parsing %v\n", problemID)
			mu.Unlock()

			err = os.MkdirAll(path, os.ModePerm)
			if err != nil {
				return
			}
			URL := fmt.Sprintf(urlFormatter, problemID)

			samples, standardIO, err := c.ParseProblem(URL, path, &mu)
			if err != nil {
				return
			}

			warns := ""
			if !standardIO {
				warns = color.YellowString("Non standard input output format.")
			}
			mu.Lock()
			if err != nil {
				color.Red("Failed %v. Error: %v", problemID, err.Error())
			} else {
				ansi.Printf("%v %v\n", color.GreenString("Parsed %v with %v samples.", problemID, samples), warns)
			}
			mu.Unlock()
		}(problemID, paths[i])
	}
	wg.Wait()
	return
}
